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

"""LongBet: time-varying heterogeneous treatment effects in panel data.

A JAX implementation of LongBet, built on ``bartz``. See :class:`LongBet`.
"""

from longbet._config import LongBetConfig
from longbet._design import Design
from longbet._diagnostics import StabilityResult, att_stability, compute_ess, compute_rhat
from longbet._encourage import (
    encouragement_effects,
    encouragement_summary,
    plot_encouragement,
    preadoption_dip_diagnostic,
    validate_encouragement,
)
from longbet._encourage_model import (
    ConditionalEffect,
    ConditionalEncouragementPrediction,
    EncouragementComparison,
    EncouragementDraws,
    EncouragementPrediction,
    IntensityConditionalPrediction,
    IntensityEncouragementPrediction,
    LongBetEncourage,
    adoption_distribution,
    offer_effect_on_outcome,
)
from longbet._encourage_design import (
    EncouragementDesign,
    EncouragementDesignResult,
    design_encouragement_effects,
)
from longbet._encourage_intensity import (
    IntensityReference,
    compose_offer_effect,
    intensity_reference,
    lag_matrix,
    validate_intensity,
)
from longbet._io import load_npz, save_npz
from longbet._multi_io import load_multi_npz, save_multi_npz
from longbet._multi_model import (
    LongBetMulti,
    LongBetMultiPrediction,
    effect_draws,
    effect_draws_from_arrays,
    joint_prob,
    outcome_correlation,
    reduce_joint_masks,
)
from longbet._rollout import (
    ROLLOUT_COLORS,
    adoption_cohorts,
    plot_rollout,
    rollout_summary,
)
from longbet._model import (
    LongBet,
    LongBetPrediction,
    available_devices,
    derive_exposure,
    get_att,
    get_catt,
    lag_dose_matrix,
    resolve_device,
)
from longbet._summary import (
    BlockAccumulator,
    PosteriorSummary,
    streaming_summary,
    summarize_draws,
)
from longbet._latecontext import LateContextSpec, late_context_evidence
from longbet._rank_k import (
    ExposureBasisSpec,
    Nested3ArmReconciliationResult,
    PlaceboAuditResult,
    RankKConfigWarning,
    RankKLongBetConfig,
    RankKLongBetResult,
    fit_rank_k_hybrid,
    fit_rank_k_longbet,
    gaussian_summary_prior_update,
    orthonormal_exposure_basis,
    prior_data_conflict_report,
    rank_k_placebo_audit,
    reconcile_nested_3arm_draws,
    reconcile_nested_3arm_monthly_draws,
)
from longbet._rank_k_plots import (
    DecayExtrapolationResult,
    PriorSensitivityResult,
    PriorTier,
    RankKPriorDraws,
    assert_nested_additivity,
    evaluate_prior_sensitivity,
    extrapolate_rank_k_decay,
    nested_3arm_trajectories,
    plot_rank_k_decay_extrapolation,
    plot_rank_k_prior_calibration,
    plot_rank_k_prior_sensitivity,
    plot_rank_k_prior_vs_posterior,
    prior_vs_posterior_summary,
    roi_break_even_table,
    sample_rank_k_prior_draws,
    standardized_prior_family_table,
    ytlr_prior_tiers,
)
from longbet._forecast import (
    CausalFeatures,
    build_causal_features,
    extract_pre_anchor_loadings,
    patch_ridge_forecaster,
)

__version__ = "0.6.0"


def build_causal_timesfm_features(*args, **kwargs):  # pragma: no cover - thin deprecated alias
    """Deprecated alias of :func:`build_causal_features` (renamed in 0.3.0)."""
    import warnings

    warnings.warn(
        "build_causal_timesfm_features was renamed to build_causal_features in longbet 0.3.0; "
        "the pretrained-TimesFM loader was removed (the built-in patch-ridge forecaster is the default)",
        DeprecationWarning,
        stacklevel=2,
    )
    return build_causal_features(*args, **kwargs)


#: Deprecated alias of :class:`CausalFeatures` (renamed in 0.3.0).
CausalTimesFMFeatures = CausalFeatures

__all__ = [
    "ROLLOUT_COLORS",
    "BlockAccumulator",
    "CausalFeatures",
    "LateContextSpec",
    "ConditionalEffect",
    "ConditionalEncouragementPrediction",
    "DecayExtrapolationResult",
    "Design",
    "EncouragementComparison",
    "EncouragementDraws",
    "EncouragementDesign",
    "EncouragementDesignResult",
    "EncouragementPrediction",
    "ExposureBasisSpec",
    "IntensityConditionalPrediction",
    "IntensityEncouragementPrediction",
    "IntensityReference",
    "LongBet",
    "LongBetConfig",
    "LongBetEncourage",
    "LongBetMulti",
    "LongBetMultiPrediction",
    "LongBetPrediction",
    "Nested3ArmReconciliationResult",
    "PosteriorSummary",
    "PriorSensitivityResult",
    "PriorTier",
    "RankKConfigWarning",
    "RankKLongBetConfig",
    "RankKLongBetResult",
    "RankKPriorDraws",
    "StabilityResult",
    "adoption_cohorts",
    "adoption_distribution",
    "assert_nested_additivity",
    "att_stability",
    "available_devices",
    "build_causal_features",
    "compose_offer_effect",
    "compute_ess",
    "compute_rhat",
    "derive_exposure",
    "design_encouragement_effects",
    "effect_draws",
    "effect_draws_from_arrays",
    "encouragement_effects",
    "encouragement_summary",
    "evaluate_prior_sensitivity",
    "extrapolate_rank_k_decay",
    "fit_rank_k_hybrid",
    "fit_rank_k_longbet",
    "gaussian_summary_prior_update",
    "get_att",
    "get_catt",
    "intensity_reference",
    "joint_prob",
    "lag_dose_matrix",
    "lag_matrix",
    "late_context_evidence",
    "load_multi_npz",
    "load_npz",
    "nested_3arm_trajectories",
    "offer_effect_on_outcome",
    "orthonormal_exposure_basis",
    "outcome_correlation",
    "patch_ridge_forecaster",
    "plot_rollout",
    "plot_encouragement",
    "plot_rank_k_decay_extrapolation",
    "plot_rank_k_prior_calibration",
    "plot_rank_k_prior_sensitivity",
    "plot_rank_k_prior_vs_posterior",
    "preadoption_dip_diagnostic",
    "prior_data_conflict_report",
    "prior_vs_posterior_summary",
    "rank_k_placebo_audit",
    "reconcile_nested_3arm_draws",
    "reconcile_nested_3arm_monthly_draws",
    "resolve_device",
    "reduce_joint_masks",
    "roi_break_even_table",
    "rollout_summary",
    "sample_rank_k_prior_draws",
    "save_multi_npz",
    "save_npz",
    "standardized_prior_family_table",
    "streaming_summary",
    "summarize_draws",
    "validate_encouragement",
    "validate_intensity",
    "ytlr_prior_tiers",
    "__version__",
]
