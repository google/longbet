# Benchmark and audit records

The Markdown files in this directory and in `encouragement/` are the historical
record of the sampler and model investigations that led to the current design:
the slow-mode diagnosis, the exact-block repairs that did not resolve it, the
variance-prior propriety proof, the tree-move validation, the parallel-tempering
measurements, and the encouragement comparison studies. They describe the
package as it was when each was written. Several of the options they mention
(adaptive coding, shared treatment trees, exposure splits, the reduced-form,
direct-smooth, orthogonal and stump-hazard encouragement engines, the collapsed
and joint-Gaussian experimental blocks) were removed in the 2026-09-14
simplification, and the scripts that exercised them were removed with them so
that nothing in the repository imports code that no longer exists.

The conclusions that shaped the current package are in
[`tempering_plan.md`](tempering_plan.md) (why the reduced-form encouragement
model was replaced by the adoption-clock model, and what tempering does and
does not buy) and [`tree_moves_review.md`](tree_moves_review.md) (the CHANGE and
REGROW moves that every sweep now runs).
