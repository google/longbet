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

"""The unified design matrix, and the spec that lets ``predict`` rebuild it.

Both forests split on **one** binned predictor matrix; which columns each may
use is controlled by its own ``max_split`` vector, with a zero entry blocking a
column.  That is how ``nu`` is denied the propensity score and the exposure
index, how ``mu`` is denied the exposure index, and how ``split_calendar_trt``
is implemented -- without a second copy of a matrix that dominates device
memory.

The layout is recorded as an ordered list of :class:`Block`s, persisted with the
model, and replayed by ``predict``.  Rebuilding from the same spec is what makes
it structurally impossible for a column present at fit time to go missing at
predict time: :func:`Design.build` fails loudly if a block's input is absent,
instead of silently producing a narrower matrix whose column indices no longer
line up with the fitted trees.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import jax
import jax.numpy as jnp
import numpy as np
from bartz.prepcovars import GivenSplitsBinner, UniqueQuantileBinner
from jaxtyping import Array, Key

#: Source key for each block, i.e. which user input it is built from.
BlockSource = str


def _binner_splits(binner: Any) -> np.ndarray:
    """Return a binner's cutpoint matrix.

    ``bartz.prepcovars.Binner`` exposes ``max_split`` and ``bin`` publicly but
    keeps the cutpoints in ``_splits``.  This is the project's only use of that
    attribute; ``tests/test_design.py`` pins the round trip through
    ``GivenSplitsBinner`` so an upstream rename is caught on the next bump.
    """
    return np.asarray(binner._splits)


@dataclass(frozen=True)
class Block:
    """One contiguous group of columns in the unified design matrix.

    Attributes
    ----------
    name
        Identifier, unique within a design.
    source
        Which user input supplies the raw values: one of ``'x'``, ``'x_trt'``,
        ``'x_tv'``, ``'x_trt_tv'``, ``'t'``, ``'s'``, ``'ps'``.
    splits
        Cutpoints, shape ``(n_cols, n_cuts)``.
    max_split
        Per-column split bound, shape ``(n_cols,)``.
    mu_visible, nu_visible
        Whether the prognostic / treatment forest may split on these columns.
    """

    name: str
    source: BlockSource
    splits: np.ndarray
    max_split: np.ndarray
    mu_visible: bool
    nu_visible: bool

    @property
    def n_cols(self) -> int:
        """Number of columns contributed by this block."""
        return int(self.max_split.shape[0])

    def bin(self, raw: np.ndarray) -> np.ndarray:
        """Bin raw values of shape ``(n_cols, n_cells)`` using this block's cuts."""
        raw = np.asarray(raw, dtype=np.float32)
        if raw.shape[0] != self.n_cols:
            raise ValueError(
                f"block {self.name!r} expects {self.n_cols} row(s), got {raw.shape[0]}"
            )
        binner = GivenSplitsBinner(
            jnp.asarray(raw), xinfo=jnp.asarray(self.splits, dtype=jnp.float32)
        )
        return np.asarray(binner.bin(jnp.asarray(raw)))

    def to_dict(self) -> dict[str, Any]:
        """Serialize to plain JSON-compatible data."""
        return {
            "name": self.name,
            "source": self.source,
            "splits": np.asarray(self.splits, dtype=np.float64).tolist(),
            "max_split": np.asarray(self.max_split, dtype=np.int64).tolist(),
            "mu_visible": bool(self.mu_visible),
            "nu_visible": bool(self.nu_visible),
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Block:
        """Rebuild from :meth:`to_dict` output."""
        return cls(
            name=d["name"],
            source=d["source"],
            splits=np.asarray(d["splits"], dtype=np.float32).reshape(
                len(d["max_split"]), -1
            ),
            max_split=np.asarray(d["max_split"], dtype=np.uint8),
            mu_visible=bool(d["mu_visible"]),
            nu_visible=bool(d["nu_visible"]),
        )


@dataclass(frozen=True)
class Design:
    """An ordered list of blocks plus the two forests' ``max_split`` vectors."""

    blocks: tuple[Block, ...]

    @property
    def max_split_mu(self) -> np.ndarray:
        """Split bounds for the prognostic forest; 0 blocks a column."""
        return np.concatenate(
            [
                b.max_split if b.mu_visible else np.zeros_like(b.max_split)
                for b in self.blocks
            ]
        ).astype(np.uint8)

    @property
    def max_split_nu(self) -> np.ndarray:
        """Split bounds for the treatment forest; 0 blocks a column."""
        return np.concatenate(
            [
                b.max_split if b.nu_visible else np.zeros_like(b.max_split)
                for b in self.blocks
            ]
        ).astype(np.uint8)

    @property
    def n_cols(self) -> int:
        """Total number of columns."""
        return sum(b.n_cols for b in self.blocks)

    def build(self, raw: dict[BlockSource, np.ndarray]) -> np.ndarray:
        """Assemble the binned design matrix, shape ``(p, n_cells)``.

        Parameters
        ----------
        raw
            Raw values per source key, each of shape ``(n_cols, n_cells)``.

        Raises
        ------
        ValueError
            If any block's source is missing. A design that was fitted with a
            propensity score must be rebuilt with one, or the fitted trees'
            column indices would no longer refer to the columns they were
            grown on.
        """
        parts = []
        for b in self.blocks:
            if b.source not in raw:
                raise ValueError(
                    f"design block {b.name!r} needs input {b.source!r}, which was not "
                    f"supplied. The model was fitted with it, so predictions require "
                    f"it too."
                )
            parts.append(b.bin(raw[b.source]))
        return np.vstack(parts)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to plain JSON-compatible data."""
        return {"blocks": [b.to_dict() for b in self.blocks]}

    def validate_panel_grids(self, n_periods: int, max_exposure: int) -> None:
        """Reject archived fits whose time predictors were silently discarded.

        In particular, early multi-outcome fits built both grids with the
        temporary scalar estimator's zero dimensions. Rebuilding cutpoints at
        prediction time cannot repair trees trained on that different design.
        """
        for name, n_cuts in (("t", n_periods - 1), ("s", max_exposure)):
            blocks = [b for b in self.blocks if b.name == name]
            expected = (np.arange(max(0, n_cuts), dtype=np.float32) + 0.5)[None, :]
            if (n_cuts < 0 or len(blocks) != 1
                    or blocks[0].source != name
                    or not np.array_equal(blocks[0].splits, expected)
                    or not np.array_equal(blocks[0].max_split, [n_cuts])):
                raise ValueError(
                    f"Invalid fitted {name!r} grid for the saved panel dimensions. "
                    "This may be a multi-outcome fit from before the time-grid "
                    "correction. Refit from the original data; loading old draws "
                    "cannot repair the model."
                )

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Design:
        """Rebuild from :meth:`to_dict` output."""
        return cls(blocks=tuple(Block.from_dict(b) for b in d["blocks"]))


def quantile_block(
    name: str,
    source: BlockSource,
    raw: np.ndarray,
    max_bins: int,
    *,
    mu_visible: bool,
    nu_visible: bool,
    key: Key[Array, ''] | None = None,
) -> Block:
    """Build a block for continuous covariates, using quantile cutpoints.

    ``UniqueQuantileBinner`` is BART's usual choice for continuous predictors.

    Above 100,000 cells it estimates the quantiles from a random subsample and
    demands a PRNG key, so one is always supplied: a panel of that size is the
    normal case here, not the exception. The resulting cutpoints are persisted
    on the block, so prediction replays them deterministically however they were
    chosen, and passing a key derived from ``random_seed`` makes the fit itself
    reproducible.
    """
    raw = np.asarray(raw, dtype=np.float32)
    if key is None:
        key = jax.random.key(0)
    binner = UniqueQuantileBinner(jnp.asarray(raw), max_bins=max_bins, key=key)
    return Block(
        name=name,
        source=source,
        splits=_binner_splits(binner),
        max_split=np.asarray(binner.max_split, dtype=np.uint8),
        mu_visible=mu_visible,
        nu_visible=nu_visible,
    )


def integer_grid_block(
    name: str,
    source: BlockSource,
    n_levels: int,
    *,
    mu_visible: bool,
    nu_visible: bool,
) -> Block:
    """Build a block for an already-discrete axis on ``0..n_levels - 1``.

    Cutpoints sit at the midpoints between consecutive integers, so every
    integer boundary is available and none is wasted. This is the right binning
    for calendar time and the exposure index; a quantile binner would place cuts
    by frequency and an even-range binner would miss boundaries.
    """
    if int(n_levels) < 1:
        raise ValueError(f"{name!r} grid requires at least one level, got {n_levels}.")
    n_cuts = int(n_levels) - 1
    splits = (np.arange(n_cuts, dtype=np.float32) + 0.5).reshape(1, -1)
    return Block(
        name=name,
        source=source,
        splits=splits,
        max_split=np.asarray([n_cuts], dtype=np.uint8),
        mu_visible=mu_visible,
        nu_visible=nu_visible,
    )
