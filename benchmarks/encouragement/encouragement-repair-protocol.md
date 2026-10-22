# Encouragement model repair: controlled comparison

This is a research protocol for the mixing failure recorded in
`benchmarks/encouragement/calibration-report.md`. It does not replace the
production model or its `calibration_status="not_established"` label.

## First experiment

Compare the existing independent Gaussian/LPM LongBet model, the same model
restricted to stumps, a Gaussian parametric direct contrast, and a direct
contrast forest with exactly sampled stump topologies. The direct models use

\[
Y_{kit}=m_k(X_i,t)+\gamma_{ki}
  +(A_i-p)1(t\ge t_0)\tau_k(X_i,t-t_0)+\epsilon_{kit}.
\]

Each leaf of the new forest contains a smooth time vector with a fixed RBF
covariance. There is no product of adaptive treatment coding, a sampled GP
trajectory, and a second treatment forest. This changes the model and prior.
The stump restriction also limits interactions; the old-model stump ablation
helps separate that restriction from the new parameterization. Exact finite
topology enumeration avoids relying on a local grow/prune proposal in this
initial experiment.

Both new models center and scale each response, use marginal baseline/effect
prior SD 1, RBF length 2 and nugget .05, independent IG(3,2) innovation and unit
variances, and independent equations. The old model receives the same variance
hyperparameters but retains its own forest/product priors: this is not an
exactly prior-matched posterior comparison. Uptake is a Gaussian working
response. It need not predict valid probabilities or absorbing adoption paths.

Use N=80, six periods, encouragement at index 2, a 40% complier share and 15%
always takers. First disable latent cross-equation unit correlation, serial
correlation and heteroskedasticity in the existing potential-outcome generator.
The richer existing `strong` scenario is a subsequent sensitivity experiment.
Potential paths precede complete balanced assignment. Two independent data seeds
and two separate sampler seeds per dataset prevent a favorable single run from
determining the conclusion. Sampler seeds do not increase the number of datasets.

Each fit uses four independent prior starts, 1,000 warmup and 2,000 retained
sweeps without thinning. Compare the first 1,000 retained sweeps with the full
2,000, preserving nesting. The same budgets in sweeps are not equal computational
budgets, so record synchronized fit time, total wall time and ESS per fit second.
The study averages encouragement contrasts over all original study units, with
the same observed covariates and weights under both assignments. Record finite
population and conditional mean truths separately, without treating horizons as
independent experiments.

## Gates and reporting

1. Verify each new conditional using independent Gaussian moment calculations,
   tiny enumerated topology probabilities, residual invariants and Geweke
   prior-data-posterior preservation checks. A deliberately broken kernel must
   fail the preservation check. Correctness precedes performance comparisons.
2. Diagnose both natural-scale ITTs and their raw Wald ratio at every reported
   horizon: rank split R-hat <=1.01, bulk and 2.5%/97.5% tail ESS >=400. Retain
   undefined ratios, negative denominators, and quantile MCSEs. Passing these
   checks is evidence about these draws, not proof of posterior convergence.
3. Report short/full ESS growth, parameter/function diagnostics and between-seed
   agreement. A minimum 1.5-fold growth is a screening flag, not a statistical
   theorem. Do not erase failures or stop selectively after a favorable seed.
4. Only advance a candidate that passes the effect checks across every initial
   data/sampler run. Investigate failed function or variance diagnostics before
   calling the sampler repaired. Coverage in two datasets cannot establish
   calibration; repeated sampler seeds must never count as independent coverage
   replications.

## Subsequent work, conditional on evidence

Advance a successful candidate to the original correlated/serial/heteroskedastic
strong scenario, then a larger independent-population calibration study. Add
joint unit intercepts and distinguish their covariance from serial innovations.
Binary/binary inference needs a separately identified multivariate probit model
with fixed latent scales, validated conditionals and natural-scale diagnostics.
An adoption model must separate baseline prevalence from hazards among units
not yet adopted, and reconstruct absorbing stocks. It does not identify an
unrestricted treatment-duration causal response from one binary instrument.

Randomization AR inversion is a separate reference improvement. Any finite
randomization exactness claim must specify the sharp homogeneous-effect null;
heterogeneous LATE guarantees require additional asymptotic conditions. A grid
of tested ratios is not a fully determined confidence set, especially its tails.

No public default, Python/R contract or fitted archive changes until the selected
model passes the appropriate sampler and calibration gates. If the first stage
is zero, a sampler cannot manufacture causal identification.

## Recorded amendments after the first experiment

The first four direct-forest runs passed every ITT/Wald diagnostic and the ESS
growth screen, while baseline and unit-intercept traces failed. This prompted
`direct_blocked`, an exact joint Gaussian update of all baseline leaf values and
unit intercepts given the topology and effect forest. Its first runs improved
baseline mixing but still left some individual intercept ESS below 400. The next
ablation, `direct_joint`, includes both forests' leaf values in that same joint
draw. All three direct variants have the same statistical posterior; old draws
and failures remain archived. Independent dense Gaussian and Geweke tests cover
each update, including deliberate broken-conditionals checks. These amendments
were triggered by the prescribed parameter diagnostics, not by interval coverage.

The joint update passed effect checks in all four initial runs, but a few unit
intercepts and aggregate tree-feature counts remained below parameter thresholds.
The retained windows were therefore extended uniformly to 4,000 and then 8,000
draws per chain, using exact saved sampler/RNG states with no additional warmup.
All shorter prefixes remain in separate result folders. At 4,000 draws, three
runs passed every parameter check and the remaining run's worst unit-intercept
R-hat was 1.0108. The final extension is a bounded check of whether further
sampling resolves that remaining disagreement, not a relaxation of the gate.
`audit_parameter_growth.py` records the parameter-level evidence; most failing
summaries had increasing ESS, while one tree-feature summary remained uneven.
