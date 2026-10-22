# Correlated-intercept sensitivity analysis, 2026-09-10 20:27 UTC

The primary protocol designated independent intercepts as the main candidate.
After inspecting its first 200 completed linear and smooth experiments, we add
a separate assessment of the repaired correlated-intercept option. This is a
secondary analysis, not a replacement or retroactive change to the primary test.

Use the same 100 datasets in each strong-stage scenario (linear, smooth, abrupt),
the same seeds, priors, tree counts, 500 warmup iterations and 1,000 saved draws
in each of four chains. Change only `correlated_intercepts=True`, retaining the
existing proper inverse-Wishart default (df 7, scale 4). Use six worker processes.
Keep every fit and the same references, accuracy summaries and pointwise coverage.
Report the primary results regardless of this analysis's outcome. Do not pool
the two configurations or choose the better one separately for each dataset.

The original accuracy, coverage and diagnostic criteria remain useful descriptive
benchmarks. Nominal paired-bootstrap intervals here are exploratory: a promising
secondary result would need independent confirmation before a superiority claim.
This analysis does not assess correlated intercepts under weak/absent relevance
or serial heteroskedastic errors and cannot support claims in those settings.
