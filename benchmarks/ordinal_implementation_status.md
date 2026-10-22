# Ordinal implementation audit

Objective: implement the full ordered-probit scope in `ordinal.md`, phases 1–5.
Status: **in progress**. Do not interpret completed fast checks as completion of
the statistical benchmark or the R suite. The original specification remains
unchanged in `ordinal.md`.

## Implementation and evidence

| Requirement | Current implementation / authoritative evidence |
|---|---|
| Declared categories, strict label/count/prior validation, missing/empty categories, correct offset | `_config.py`, `_ordinal.prepare_ordinal`; `test_ordinal_primitives.py` |
| Proper ordered-normal prior; sequential Gibbs and marginalized partially collapsed MH | `_ordinal.sample_cutpoints`, `_ordinal.sample_cutpoints_marginalized`; SciPy conditional and ordered-prior checks |
| Stable float32 interval draws, jit/vmap, extreme/narrow bounds and explicit failures | `_ordinal.sample_truncated_normal`; independent SciPy probability-integral-transform tests. Stock JAX 0.11.1 degenerated at [8,9] and [12,infinity], so it is not used for K>2. |
| Conditional latent mean/scale, missing placeholders, residual invariants, unit variance | `_state.py`, `_step.py`; scalar invariant tests and SUR conditional-normal oracle |
| Per-chain cutpoints and shared inputs; legacy RNG schedule; K=2 equality | Fold-in tags 8101–8103 for initialization and 8111–8112 for sweeps; legacy split counts unchanged; exact binary comparisons include absent classes, missingness, multiple chains |
| Aligned threshold traces without per-sweep latents | `_loop.py`, `_multi_loop.py`; saved-sweep/thinning and multi-shape tests |
| Paired category probabilities/ATT/custom scores with cells first, draws last | `_model.py`; algebra, aggregation, score, shape and unsupported-interface tests |
| Exact summaries with a panel-independent working memory budget | `_summary.py` and the prediction cell loop; allocation/transform instrumentation with cache disabled |
| Scalar NPZ schema and legacy/stale archive handling | `_io.py`; round-trip, corruption and old nonordinal archive tests |
| Heterogeneous K in user/internal order, discrete-first stable partition | `_multi_input.py`, `_multi_state.py`, `_multi_model.py`; mixed continuous/ordinal/binary/ordinal tests |
| Identified triangular SUR, downstream precision, missingness restrictions, shared trees | Existing `_sur.py` and `_shared_forest.py` reused with ordinal child sweeps; conditional distribution, residual, topology and independent-scalar oracles |
| Multi format 3, parent and extracted child persistence | `_multi_io.py`; shared and unshared round-trips and malformed order/count/cutpoint/variance/loading tests |
| R counts, labels, probability methods, category axes, rehydration | `R/ordinal.R` and scalar/multi wrappers; `testthat/test-ordinal.R`. Full R CMD check passed with 385 assertions, no failures/warnings/skips, including actual fresh-process rehydration. |
| Public docs and shared API contract | README ordinal section, `man/*.Rd`, `NAMESPACE`, both identical YAML contracts; Python contract checks pass |
| Required slow panel, ten-seed coverage, held-out calibration, intercept/missing extension, diagnostics/runtime | `ordinal_validation.py`, `report_ordinal_validation.py`, `test_ordinal_validation.py`; first full benchmark completed and inspected. Remaining seeds and final report pending. |

For K=2, exact binary compatibility includes the existing binary kernel's unused
random missing-cell latents. The K>2 path zeros missing latents. Both exclude
missing cells from likelihoods, residuals and threshold bounds.

## Completed commands

These results were read from completed process output, not inferred from the
existence of log files:

- `pytest tests/test_ordinal_primitives.py -q`: **55 passed**.
- `pytest tests/test_step.py tests/test_chain_axis.py tests/test_model.py tests/test_ordinal_primitives.py -q`: **92 passed**, 160.97 seconds.
- `pytest tests/test_ordinal_model.py -q`: **18 passed**, 131.29 seconds. This included the temporary explicit joint rejection test, removed when joint implementation was added.
- `pytest tests/test_ordinal_predict.py -q`: **32 passed**, 47.70 seconds.
- `pytest tests/test_ordinal_multi.py -xq`: **41 passed**, 125.50 seconds.
- `pytest tests/test_multi_inputs.py tests/test_multi_scalar_equivalence.py tests/test_multi_io.py tests/test_multi_saving.py tests/test_multi_loading.py -q`: **31 passed**, 7,741.51 seconds, one expected always-treated-unit warning.
- `pytest tests/test_contract.py -q`: **7 passed**, 3.82 seconds.
- `pytest tests/test_ordinal_validation.py -m 'not slow' -q`: **1 passed**, 11 deselected, 4.38 seconds.
- Final worktree: `pytest tests/test_ordinal_model.py tests/test_ordinal_predict.py tests/test_ordinal_multi.py tests/test_ordinal_primitives.py tests/test_multi_inputs.py tests/test_contract.py -q`: **181 passed**, 389.83 seconds. This includes the public all-top-category fit, scalar K=2 archive corruption, explicit paired forest-input/offset tests, and continuous missing-value sanitization.
- `pytest tests/test_ordinal_multi.py -k two_category_joint_archives -q`: **2 passed**, 41 deselected, 39.95 seconds. Both shared and unshared K=2 multi archives preserve `(draws,0)` thresholds and reject an absent empty array.
- `R CMD build` followed by `R CMD check --no-manual`: **Status: OK**, all documentation and example checks passed; **385 assertions, 0 failures, 0 R warnings, 0 skips**, 347.10 seconds for the test runner. This includes the whole ordinal suite, full/summary Python/R parity, singleton category axes, and fresh-process parent/child/prediction rehydration. The expected Python singleton-unit warning appears on stderr.

Python commands use `.venv/bin/python -m pytest`. The host `/tmp` reached its
per-user quota while other work was present. New scratch files and logs use
`benchmarks/ordinal_results/.tmp` and `benchmarks/ordinal_results/logs`, both
ignored by git. R tests use the existing `longbet-rtest` Docker image with this
worktree mounted at `/src` and `PYTHONPATH=/src/src`. CPU memory contention has
caused long pauses; benchmark timings retain those delays. R CMD build copies
before filtering, so package staging applied the existing `.Rbuildignore` rules
first to avoid copying excluded benchmark data and host-absolute pytest symlinks
into the container. The checked package used the current R/man/inst/tests files.
R versions: R 4.6.1, reticulate 1.47.0, testthat 3.3.2; its Python engine used
Python 3.12.3, JAX/jaxlib 0.11.1, bartz 0.12.1, equinox 0.13.8, NumPy 2.5.3,
SciPy 1.18.1, and ArviZ 1.3.0. The statistical benchmark's NumPy version is
2.5.2; all of its versions are saved per dataset.

## Remaining gates

1. Complete seeds 20260910–20260918 and the intercept/missing
   panel. Generate the coverage report with binomial uncertainty and run all
   eleven slow artifact/validation tests. The partially collapsed marginalized
   threshold proposal resolved the cutpoint mixing bottleneck (cutpoint max R-hat
   dropped from 3.11 to 1.01, bulk ESS improved from 4.5 to 347). Category and score
   effects converge smoothly across chains.
2. Inspect all final artifacts and diffs against every section of `ordinal.md`;
   update the report and mark the goal complete only when all required evidence
   is present. The final 183 Python checks and the full R package check passed.

Section 10 defines a separate cloglog evaluation before replacing ordered
probit. This implementation retains probit and exposes no unused link option.
No claim of cloglog superiority, conjugacy for GP-weighted treatment leaves,
or production support is made.
