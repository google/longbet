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

"""Posterior summaries, including a streaming path for large panels.

An ``N x T x draws`` array is 12 GB in float32 at 100,000 x 30 x 1000.
:func:`streaming_summary` therefore evaluates the posterior in blocks of *cells*
and reduces each block to its summary immediately, so peak memory is
``block_size * draws`` rather than ``N * T * draws``.  Blocking over cells rather
than over draws keeps the quantiles **exact** -- every draw for a given cell is
present when that cell is reduced.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import NamedTuple

import numpy as np

#: Target elements held in memory per block; ~32 MB in float32.
_TARGET_BLOCK_ELEMENTS = 8_000_000


def choose_ordinal_block_size(n_draws: int, n_cells: int, num_categories: int) -> int:
    """Bound ordinal working storage to ~32 MB, including float64 CDF buffers.

    Budget sixteen float64 category buffers plus sixteen float32 latent/forest
    buffers per cell/draw. This includes paired arm probabilities, CDF bounds,
    log-CDFs, quantile work, and their transient copies. Required retained
    summaries and exposure ATT draws are outside this working budget.
    """
    per_cell = max(1, n_draws) * (16 * 8 * num_categories + 16 * 4)
    return min(max(1, n_cells), max(1, _TARGET_BLOCK_ELEMENTS * 4 // per_cell))


class PosteriorSummary(NamedTuple):
    """Posterior point estimates and credible interval bounds."""

    mean: np.ndarray
    std: np.ndarray
    lower: np.ndarray
    upper: np.ndarray


def summarize_draws(
    draws: np.ndarray,
    alpha: float = 0.05,
    axis: int = 0,
) -> PosteriorSummary:
    """Summarize an in-memory array of draws along ``axis``.

    Suitable for small quantities (the ATT series, ``beta``). For per-cell
    quantities on a large panel use :func:`streaming_summary`.
    """
    arr = np.asarray(draws)
    q_low = 100.0 * (alpha / 2.0)
    q_high = 100.0 * (1.0 - alpha / 2.0)
    lower, upper = np.percentile(arr, [q_low, q_high], axis=axis)
    return PosteriorSummary(
        mean=np.mean(arr, axis=axis),
        std=np.std(arr, axis=axis),
        lower=lower,
        upper=upper,
    )


def choose_block_size(n_draws: int, n_cells: int) -> int:
    """Cells per block so that a block holds about 32 MB in float32."""
    if n_draws <= 0:
        return max(1, n_cells)
    return int(min(max(1, _TARGET_BLOCK_ELEMENTS // n_draws), max(1, n_cells)))


class BlockAccumulator:
    """Accumulate exact per-cell summaries one block of cells at a time.

    Blocking over *cells* rather than over draws is what keeps the quantiles
    exact: every draw for a given cell is present when that cell is reduced.
    Peak memory is ``block_size * draws`` instead of ``n_cells * draws``.
    """

    def __init__(self, n_cells: int, alpha: float = 0.05, *, dtype=np.float32) -> None:
        self.mean = np.empty(n_cells, dtype=dtype)
        self.std = np.empty(n_cells, dtype=dtype)
        self.lower = np.empty(n_cells, dtype=dtype)
        self.upper = np.empty(n_cells, dtype=dtype)
        self._q_low = 100.0 * (alpha / 2.0)
        self._q_high = 100.0 * (1.0 - alpha / 2.0)

    def update(self, lo: int, hi: int, block: np.ndarray) -> None:
        """Reduce one ``(draws, cells)`` block into the running outputs."""
        self.mean[lo:hi] = block.mean(axis=0)
        self.std[lo:hi] = block.std(axis=0)
        low, high = np.percentile(block, [self._q_low, self._q_high], axis=0)
        self.lower[lo:hi] = low
        self.upper[lo:hi] = high

    def result(self) -> PosteriorSummary:
        """Return the accumulated summary."""
        return PosteriorSummary(self.mean, self.std, self.lower, self.upper)


def streaming_summary(
    block_fn: Callable[[int, int], np.ndarray],
    n_cells: int,
    n_draws: int,
    alpha: float = 0.05,
    block_size: int | None = None,
) -> PosteriorSummary:
    """Summarize a per-cell posterior without materializing all draws.

    Parameters
    ----------
    block_fn
        ``block_fn(lo, hi)`` returns the draws for cells ``[lo, hi)`` as an
        array of shape ``(n_draws, hi - lo)``.
    n_cells
        Total number of cells to summarize.
    n_draws
        Number of posterior draws per cell.
    alpha
        Tail probability; the interval is ``[alpha/2, 1 - alpha/2]``.
    block_size
        Cells per block. Defaults to :func:`choose_block_size`.

    Returns
    -------
    PosteriorSummary
        Flat arrays of length ``n_cells``. Quantiles are exact.
    """
    if block_size is None:
        block_size = choose_block_size(n_draws, n_cells)
    block_size = max(1, int(block_size))

    acc = BlockAccumulator(n_cells, alpha)
    for lo in range(0, n_cells, block_size):
        hi = min(lo + block_size, n_cells)
        acc.update(lo, hi, np.asarray(block_fn(lo, hi)))
    return acc.result()
