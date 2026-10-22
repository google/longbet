# Copyright 2026 Google LLC

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     https://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""The Python side of the two-front-door contract (contract/longbet-api.yaml).

Checked in **both** directions: every contract entry must match the dataclass,
and every dataclass field must appear in the contract or be declared internal.
A one-directional check cannot catch a new option that was never written down.
"""

import dataclasses
from pathlib import Path

import pytest
import yaml

from longbet import LongBetConfig

CONTRACT_PATH = Path(__file__).resolve().parents[1] / "contract" / "longbet-api.yaml"


@pytest.fixture(scope="module")
def contract():
    assert CONTRACT_PATH.exists(), f"contract file not found: {CONTRACT_PATH}"
    with open(CONTRACT_PATH, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _defaults():
    """LongBetConfig defaults as the plain model resolves them: a model-dependent
    option is None in the dataclass and MODEL_DEFAULTS for longbet()."""
    from longbet._config import MODEL_DEFAULTS
    return {f.name: MODEL_DEFAULTS[f.name] if f.default is None and f.name in MODEL_DEFAULTS
            else f.default for f in dataclasses.fields(LongBetConfig)}


def test_model_dependent_options_match_the_config(contract):
    from longbet._config import MODEL_DEFAULTS
    flagged = {a["python_name"] for a in contract["arguments"] if a.get("model_dependent")}
    assert flagged == set(MODEL_DEFAULTS)
    fields = {f.name: f.default for f in dataclasses.fields(LongBetConfig)}
    assert all(fields[name] is None for name in flagged)


def test_encouragement_components_match_the_engine(contract):
    from longbet._encourage_model import DEFAULT_N_SKIP, MODEL_COMPONENTS
    api = contract["encouragement_model_api"]
    assert {k for k, v in api["components"].items() if v is True} == set(MODEL_COMPONENTS)
    assert api["defaults"] == {"n_skip": DEFAULT_N_SKIP}
    assert set(api["removed_options"]) == {"linear_trt_horseshoe", "linear_trt_outcome_only",
                                           "adoption_frailty"}


def test_every_contract_argument_matches_the_config(contract):
    defaults = _defaults()
    for arg in contract["arguments"]:
        name = arg["python_name"]
        assert name in defaults, f"contract lists {name!r}, which LongBetConfig lacks"
        assert defaults[name] == arg["default"], (
            f"default mismatch for {name}: contract says {arg['default']!r}, "
            f"LongBetConfig says {defaults[name]!r}"
        )


def test_every_config_field_is_in_the_contract(contract):
    documented = {a["python_name"] for a in contract["arguments"]}
    internal = set(contract.get("internal", []))
    undocumented = set(_defaults()) - documented - internal
    assert not undocumented, (
        f"LongBetConfig fields missing from the contract: {sorted(undocumented)}. "
        "Add them to contract/longbet-api.yaml (and to the R front door), or list "
        "them under `internal`."
    )


def test_r_names_are_unique(contract):
    r_names = [a["r_name"] for a in contract["arguments"] if a.get("expose_in_r")]
    assert len(r_names) == len(set(r_names)), "duplicate r_name in the contract"


def test_yaml_files_are_byte_identical():
    c1 = (Path(__file__).resolve().parents[1] / "contract" / "longbet-api.yaml").read_bytes()
    c2 = (Path(__file__).resolve().parents[1] / "inst" / "contract" / "longbet-api.yaml").read_bytes()
    assert c1 == c2, "contract/longbet-api.yaml and inst/contract/longbet-api.yaml differ"


def test_multi_contract_structure(contract):
    assert "multi_api" in contract, "contract missing multi_api section"
    m_api = contract["multi_api"]
    from longbet._sur import SAMPLER_SEMANTICS
    from longbet._io import PRECISION_CACHE_VERSION
    assert m_api.get("sampler_semantics") == SAMPLER_SEMANTICS
    assert m_api.get("precision_cache_version") == PRECISION_CACHE_VERSION
    assert "proper by construction" in m_api["variance_prior_requirement"]
    assert "inputs" in m_api
    assert "fit_returns" in m_api
    assert "predict_returns" in m_api
    assert "utilities" in m_api
    for util in ("effect_draws", "joint_prob", "outcome_correlation"):
        assert util in m_api["utilities"]


def test_python_multi_exports():
    import longbet
    for symbol in (
        "LongBetMulti",
        "LongBetMultiPrediction",
        "effect_draws",
        "joint_prob",
        "outcome_correlation",
        "effect_draws_from_arrays",
        "reduce_joint_masks",
    ):
        assert hasattr(longbet, symbol), f"longbet missing export {symbol}"
        assert symbol in longbet.__all__, f"{symbol} missing from longbet.__all__"


def test_encouragement_contract_signatures_exports_and_columns(contract):
    import inspect
    import numpy as np
    import longbet

    api = contract["encouragement_api"]
    namespace = (CONTRACT_PATH.parents[1] / "NAMESPACE").read_text()
    for name, spec in api["functions"].items():
        assert name in longbet.__all__
        assert f"export({name})" in namespace
        signature = inspect.signature(getattr(longbet, name))
        expected = set(spec["required"]) | set(spec["defaults"]) | set(spec.get("python_only_defaults", {}))
        assert set(signature.parameters) == expected
        for key in spec["required"]:
            assert signature.parameters[key].default is inspect.Parameter.empty
        for key, value in {**spec["defaults"], **spec.get("python_only_defaults", {})}.items():
            assert signature.parameters[key].default == value

    z = np.array([[0, 1], [0, 1], [0, 0], [0, 0]])
    assert list(longbet.validate_encouragement(z, z)) == api["metadata_fields"]
    assert list(longbet.encouragement_summary(z, z)) == api["summary_columns"]
    result = longbet.encouragement_effects(z, z, z)
    assert list(result) == api["effects_columns"]
    assert result.inference.iloc[0] == api["inference"]["method"]
    assert result.wald_set_type.iloc[0] in api["inference"]["set_types"]


def test_r_frontdoor_signatures_and_defaults(contract):
    """Ensure R functions longbet and longbet_multi have formals matching the contract."""
    import re
    repo_root = Path(__file__).resolve().parents[1]

    # Parse R/longbet.R
    lb_r = (repo_root / "R" / "longbet.R").read_text(encoding="utf-8")
    lb_sig_match = re.search(r"longbet\s*<-\s*function\s*\((.*?)\)\s*\{", lb_r, re.DOTALL)
    assert lb_sig_match, "Could not find longbet signature in R/longbet.R"
    lb_params = [p.strip().split("=")[0].strip() for p in lb_sig_match.group(1).split(",") if p.strip()]

    # Parse R/multi.R
    multi_r = (repo_root / "R" / "multi.R").read_text(encoding="utf-8")
    multi_sig_match = re.search(r"longbet_multi\s*<-\s*function\s*\((.*?)\)\s*\{", multi_r, re.DOTALL)
    assert multi_sig_match, "Could not find longbet_multi signature in R/multi.R"
    multi_params = [p.strip().split("=")[0].strip() for p in multi_sig_match.group(1).split(",") if p.strip()]

    # Check exposed arguments in both functions
    for arg in contract["arguments"]:
        if not arg.get("expose_in_r"):
            continue
        r_name = arg["r_name"]
        assert r_name in lb_params, f"longbet in R/longbet.R missing formal '{r_name}'"
        assert r_name in multi_params, f"longbet_multi in R/multi.R missing formal '{r_name}'"

    # Check NAMESPACE exports
    ns = (repo_root / "NAMESPACE").read_text(encoding="utf-8")
    for exp in (
        "export(longbet_multi)",
        "export(effect_draws)",
        "export(joint_prob)",
        "export(outcome_correlation)",
        "S3method(predict, longbet_multi)",
        "S3method(print, longbet_multi)",
        "S3method(summary, longbet_multi)",
        "S3method(print, longbet_multi.pred)",
        'S3method("[", longbet_multi)',
        'S3method("[", longbet_multi.pred)',
    ):
        assert exp in ns, f"NAMESPACE missing {exp}"


def test_encouragement_model_contract(contract):
    import inspect
    import longbet
    from longbet._encourage_model import ARCHIVE_VERSION, INFERENCE_VERSION

    api = contract["encouragement_model_api"]
    assert api["archive_version"] == ARCHIVE_VERSION
    assert api["inference_version"] == INFERENCE_VERSION
    assert api["calibration_status"] == "not_established"
    cls = longbet.LongBetEncourage
    for name in ("LongBetEncourage", "EncouragementPrediction", "EncouragementComparison"):
        assert name in longbet.__all__
    init = inspect.signature(cls)
    for key in api["constructor"]["required"]:
        assert init.parameters[key].default is inspect.Parameter.empty
    for key, value in api["constructor"]["defaults"].items():
        assert init.parameters[key].default == value
    model = cls()
    for key, value in api["constructor"]["wrapper_default_priors"].items():
        assert getattr(model.config, key) == value
    for method in ("fit", "predict", "bootstrap_comparison"):
        sig = inspect.signature(getattr(cls, method))
        spec = api[method]
        expected = {"self"} | set(spec.get("required", [])) | set(spec["defaults"])
        assert set(sig.parameters) == expected
        for key in spec.get("required", []):
            assert sig.parameters[key].default is inspect.Parameter.empty
        for key, value in spec["defaults"].items():
            assert sig.parameters[key].default == value
    for method, defaults in api["posterior_methods"].items():
        sig = inspect.signature(getattr(longbet.EncouragementPrediction, method))
        assert set(sig.parameters) == {"self"} | set(defaults)
        for key, value in defaults.items():
            assert sig.parameters[key].default == value


def test_encouragement_design_contract(contract):
    import inspect
    import longbet
    api = contract["encouragement_design_api"]
    fields = {f.name: f.default for f in dataclasses.fields(longbet.EncouragementDesign)}
    assert fields == api["defaults"]
    sig = inspect.signature(longbet.design_encouragement_effects)
    assert set(sig.parameters) == set(api["required"]) | set(api["effects_defaults"])
    for key, value in api["effects_defaults"].items():
        assert sig.parameters[key].default == value
    for key in api["required"]:
        assert sig.parameters[key].default is inspect.Parameter.empty
    for name in ("EncouragementDesign", "EncouragementDesignResult", "design_encouragement_effects"):
        assert name in longbet.__all__


def test_encouragement_intensity_contract(contract):
    import inspect
    import longbet
    from longbet import _encourage_intensity as ei
    from longbet._encourage_model import (
        ARCHIVE_VERSION_INTENSITY, DEFAULT_IDENTIFICATION, ENGINE_INTENSITY, INFERENCE_VERSION_INTENSITY,
    )

    api = contract["encouragement_intensity_api"]
    assert api["engine"] == ENGINE_INTENSITY
    assert api["inference_version"] == INFERENCE_VERSION_INTENSITY
    assert api["archive_version"] == ARCHIVE_VERSION_INTENSITY
    assert tuple(api["identification_tiers"]) == ei.IDENTIFICATION_TIERS
    assert api["default_identification"] == DEFAULT_IDENTIFICATION
    assert api["tier_labels"] == ei.TIER_LABELS
    assert tuple(api["stage2_backends"]) == ei.STAGE2_BACKENDS
    assert api["calibration_status"] == "not_established"
    namespace = (CONTRACT_PATH.parents[1] / "NAMESPACE").read_text()
    for name, spec in api["reference_functions"].items():
        assert name in longbet.__all__
        assert f"export({name})" in namespace
        signature = inspect.signature(getattr(longbet, name))
        assert set(signature.parameters) == set(spec["required"]) | set(spec["defaults"])
        for key in spec["required"]:
            assert signature.parameters[key].default is inspect.Parameter.empty
        for key, value in spec["defaults"].items():
            assert signature.parameters[key].default == value
    for name in ("IntensityEncouragementPrediction", "IntensityConditionalPrediction", "IntensityReference",
                 "lag_matrix", "compose_offer_effect"):
        assert name in longbet.__all__
    model_api = contract["encouragement_model_api"]
    assert model_api["constructor"]["defaults"]["treatment"] == "absorbing"
    assert model_api["constructor"]["defaults"]["identification"] is None
    assert model_api["fit"]["defaults"]["strata"] is None
    assert set(model_api["r_fit"]["defaults"]) >= {"treatment", "identification", "lag_horizon", "demean_pre",
                                                   "stage2_backend", "strata"}
