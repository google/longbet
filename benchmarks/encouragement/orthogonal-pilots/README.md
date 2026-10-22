# Prediction-only sampling-budget pilots

These three disjoint-seed pilots compare two LongBet sampling budgets on held-out
prediction, before any comparative IV evaluation. They do not score IV truths or
competing learners. The selected budget is two chains, 200 warmup and 400 retained
iterations per chain, with the same three-candidate inner validation procedure.

The exact source used by the pilots is archived in `_iv_nuisance_pilot_source.py`.
Its hash is recorded in every JSON. The evaluated library differs only in the
default budget values; both budgets were passed explicitly in the pilot runner.
The runner is preserved as executed, including local source and output paths.
`pilot-summary.json` records the selection rationale and observed limitations.

Sampling agreement is a numerical prediction diagnostic, not a proof of MCMC
convergence or valid causal coverage. In particular, adoption predictions in the
abrupt pilot retained noticeable Monte Carlo variation.
