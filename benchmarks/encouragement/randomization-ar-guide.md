# Randomization AR research reference

`randomization_ar.py` evaluates studentized randomization tests for a supplied
grid of candidate encouragement-Wald ratios. It supports complete individual
randomization, one encouragement wave, a permanent holdout, complete outcomes
and absorbing binary adoption. The observed assignment counts define the exact
randomization law. Clustered, blocked and staggered assignment require their own
randomization mechanisms and are rejected by this research interface.

At each candidate ratio, the method tests the transformed outcome `Y - beta*D`.
Small designs are enumerated; larger designs use reproducible independent
random assignments and the conservative plus-one Monte Carlo p-value. It never
divides by the estimated first stage. The same random assignments are reused
across the grid and horizons. P-values are pointwise, not simultaneous.

Finite-sample validity requires the horizon-specific **sharp null** that every
unit's transformed potential outcome is unchanged by assignment. Exclusion and
a constant effect of current adoption imply that null. Heterogeneous effects or
treatment-history effects generally do not. The studentized heterogeneous-LATE
result in [Aronow, Chang and Lopatto](https://arxiv.org/abs/2404.18786) is
asymptotic under the paper's assumptions; this implementation makes no finite
exactness claim for heterogeneous LATE.

Use from this directory (or add it to Python's module search path):

```python
import numpy as np
from randomization_ar import randomization_ar

result = randomization_ar(y, d, z, t, beta=np.linspace(-3, 5, 161),
                          method="monte_carlo", permutations=9999, seed=42)
print(result.table)
print(result.metadata)
```

The table describes **evaluated grid points only**. Accepted endpoints do not
establish unboundedness; acceptance between grid points and in both tails remains
unknown. The `mc_p_lower/upper` columns are 95% Clopper–Pearson bounds on the ideal
enumeration tail probability, describing simulation error rather than causal
parameter uncertainty. Constant transformed outcomes produce p=1; positive
contrasts with zero estimated variance produce infinite statistics with inclusive
tie handling.

Tests enumerate an entire eight-unit randomized experiment under a sharp null
and verify conservative rejection at several nominal levels. They also cover
Monte Carlo reproducibility, the plus-one formula, infinite ties, zero first
stages, affine invariance, and unsupported requests. This is a research reference,
not a replacement for the package's documented analytic Fieller/AR API.
