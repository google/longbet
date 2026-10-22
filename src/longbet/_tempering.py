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

"""Parallel tempering and annealed burn-in for the scalar LongBet sampler.

The tempered target at inverse temperature ``beta`` is ``p(theta) L(theta)^beta``.
For the Gaussian outcome ``L^beta`` is a Gaussian likelihood with precision
``beta / sigma2`` up to a factor that depends only on ``sigma2``, so every mean,
latent and coefficient update of the ordinary sweep is tempered exactly by
running it with a per-cell conditional precision of ``beta / sigma2`` (the path
the SUR coupling already uses), and ``sigma2`` is then drawn from
``IG(a + beta m / 2, b + beta SSR / 2)``. For binary and ordinal outcomes the
latent augmentation becomes ``N(z | f, 1 / beta) 1{sign}``, a valid bridging
family whose ``beta = 1`` member is the posterior.

Replica ``c`` of a ``K``-level ladder has level ``c % K``; level 0 is the
posterior and the levels are equally spaced in ``sqrt(beta)`` down to
``beta_min``. After
every sweep adjacent levels of each ladder attempt an exchange (even pairs on
even sweeps, odd pairs on odd sweeps) accepted with probability
``min(1, exp((beta_i - beta_j)(E_j - E_i)))`` where ``E`` is the untempered
log-likelihood of the replica's state. The states move, the temperatures stay.
"""

from __future__ import annotations

import warnings
from typing import Any

import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
from jax import random
from jaxtyping import Array, Bool, Float32, Int32, Key

from longbet._state import LongBetState, split_chain_fields


def ladder_temperatures(num_replicas: int, levels: int, beta_min: float) -> Float32[Array, ' num_replicas']:
    """Inverse temperatures for ``num_replicas`` replicas laid out ladder-major."""
    # Uniform in sqrt(beta): the log-likelihood's posterior spread grows like
    # 1/sqrt(beta) as the likelihood flattens, so equal steps in sqrt(beta)
    # give roughly equal exchange rates at every level (measured on the chapter
    # panel: geometric spacing accepted 1-3 percent near beta = 1 and 40
    # percent at the hot end; sqrt spacing is even).
    root = np.linspace(1.0, np.sqrt(float(beta_min)), int(levels)) if levels > 1 else np.ones(1)
    ladder = root ** 2
    level = np.arange(num_replicas) % int(levels)
    return jnp.asarray(ladder[level], dtype=jnp.float32)


def chain_energy(state: LongBetState) -> Float32[Array, ' chains']:
    """Untempered log-likelihood of every replica, up to constants shared by all.

    Continuous: ``-(m/2) log sigma2 - SSR / (2 sigma2)``. Binary and ordinal:
    ``-SSR / 2`` for the augmented latent response (unit variance).
    """
    resid = jnp.where(state.obs_mask, state.resid, 0.0)
    if state.ar1:
        # AR(1) innovations: -(m/2) log sigma2 + (N/2) log(1 - phi^2) - |T R|^2 / (2 sigma2).
        from longbet._ar1 import transform
        innov = jax.vmap(transform, in_axes=(0, 0, None))(resid, state.ar1_phi, state.ar1_first)
        ssr = jnp.sum(jnp.square(innov), axis=-1)
        m_obs = jnp.sum(state.obs_mask.astype(jnp.float32))
        n_units = jnp.sum(state.ar1_first.astype(jnp.float32))
        return (-0.5 * m_obs * jnp.log(state.sigma2)
                + 0.5 * n_units * jnp.log1p(-jnp.square(state.ar1_phi))
                - ssr / (2.0 * state.sigma2))
    ssr = jnp.sum(jnp.square(resid), axis=-1)
    if state.num_categories >= 2:
        return -0.5 * ssr
    m_obs = jnp.sum(state.obs_mask.astype(jnp.float32))
    return -0.5 * m_obs * jnp.log(state.sigma2) - ssr / (2.0 * state.sigma2)


def swap_step(
    key: Key[Array, ''], state: LongBetState, sweep_index: Int32[Array, '']
) -> tuple[LongBetState, Bool[Array, ' chains']]:
    """One round of adjacent-level exchanges within every ladder.

    Returns the new state and, for every replica, whether it was the lower
    level of an accepted exchange with the next level.
    """
    levels = state.tempering_levels
    C = state.num_chains
    c = jnp.arange(C)
    level = c % levels
    E = chain_energy(state)
    beta = state.temperature
    parity = sweep_index % 2
    is_lead = (level % 2 == parity) & (level < levels - 1)
    nxt = jnp.minimum(c + 1, C - 1)
    log_ratio = (beta - beta[nxt]) * (E[nxt] - E)
    u = random.uniform(key, (C,), dtype=jnp.float32)
    accept_lead = is_lead & (jnp.log(u) < log_ratio)
    accept_follow = jnp.concatenate([jnp.zeros(1, bool), accept_lead[:-1]])
    perm = jnp.where(accept_lead, nxt, c)
    perm = jnp.where(accept_follow, c - 1, perm)
    per, shared = split_chain_fields(state)
    per = jax.tree_util.tree_map(lambda a: a[perm], per)
    new_state = eqx.combine(per, shared)
    # The states moved; the temperatures belong to the levels.
    new_state = eqx.tree_at(lambda s: s.temperature, new_state, beta)
    return new_state, accept_lead


#: Exchange acceptance below which a fit warns that tempering is not mixing
#: the ladder (the ladder is designed for a rate near a third).
SWAP_RATE_FLOOR = 0.05


def swap_rates(swap_accepts: Any, num_chains: int, levels: int, total_sweeps: int) -> np.ndarray:
    """Exchange acceptance per chain and adjacent level pair, ``(num_chains,
    levels - 1)``, from the per-replica counts of :func:`swap_step`."""
    # Level l leads an exchange on the sweeps whose parity matches l's.
    attempts = np.array([(total_sweeps + 1) // 2 if lvl % 2 == 0 else total_sweeps // 2
                         for lvl in range(levels - 1)], dtype=np.float64)
    accepts = np.asarray(swap_accepts, dtype=np.float64).reshape(num_chains, levels)[:, :-1]
    return accepts / np.maximum(attempts, 1.0)


def warn_if_not_exchanging(rates: np.ndarray, stacklevel: int = 3) -> None:
    """Warn when some adjacent pair of a ladder almost never exchanges."""
    if rates.size == 0 or rates.min() >= SWAP_RATE_FLOOR:
        return
    chain, pair = np.unravel_index(np.argmin(rates), rates.shape)
    warnings.warn(
        f"parallel-tempering replicas rarely exchange: levels {pair} and {pair + 1} of "
        f"chain {chain} swapped in {rates[chain, pair]:.1%} of attempts (the lowest pair; "
        f"see swap_rates_). The hot replicas then do not reach the cold chain and "
        f"tempering only adds cost: add levels (tempering_levels), raise "
        f"tempering_beta_min, or fit without tempering.", UserWarning, stacklevel=stacklevel)


def select_replicas(state: LongBetState, idx: Any) -> LongBetState:
    """Restrict every chain-axis field to the replicas ``idx``."""
    per, shared = split_chain_fields(state)
    per = jax.tree_util.tree_map(lambda a: a[idx], per)
    return eqx.combine(per, shared)


# ---------------------------------------------------------------------------
# Coupled multi-outcome sampler
# ---------------------------------------------------------------------------

def multi_chain_energy(multi_state: Any) -> Float32[Array, ' chains']:
    """Untempered SUR log-likelihood of every replica, up to shared constants.

    With unit-lower-triangular ``B = I - Gamma`` the raw residual vector has
    density ``prod_m N((B r)_m | 0, sigma2_m)`` (binary and ordinal equations
    have unit variance and their augmented latent response).
    """
    states = multi_state.states
    raw = jnp.stack([st.resid for st in states], axis=-2)          # (*C, M, n)
    Gamma = multi_state.gamma_loadings                              # (*C, M, M)
    innov = raw - jnp.einsum('...ij,...jn->...in', Gamma, raw)
    energy = 0.0
    for m, st in enumerate(states):
        e = jnp.where(st.obs_mask, innov[..., m, :], 0.0)
        ssr = jnp.sum(jnp.square(e), axis=-1)
        if st.num_categories >= 2:
            energy = energy - 0.5 * ssr
        else:
            m_obs = jnp.sum(st.obs_mask.astype(jnp.float32))
            energy = energy - 0.5 * m_obs * jnp.log(st.sigma2) - ssr / (2.0 * st.sigma2)
    return energy


def multi_set_temperature(multi_state: Any, temperature: Any) -> Any:
    """Give every child equation the same per-replica inverse temperatures."""
    temps = tuple(jnp.asarray(temperature, dtype=jnp.float32) for _ in multi_state.states)
    return eqx.tree_at(lambda s: tuple(st.temperature for st in s.states), multi_state, temps)


def multi_select_replicas(multi_state: Any, idx: Any) -> Any:
    """Restrict every chain-axis field of a coupled state to the replicas ``idx``."""
    from longbet._multi_state import split_multi_chain_fields
    per, shared = split_multi_chain_fields(multi_state)
    per = jax.tree_util.tree_map(lambda a: a[idx], per)
    return eqx.combine(per, shared)


def multi_swap_step(
    key: Key[Array, ''], multi_state: Any, sweep_index: Int32[Array, '']
) -> tuple[Any, Bool[Array, ' chains']]:
    """Adjacent-level exchanges of whole coupled states; see :func:`swap_step`."""
    from longbet._multi_state import split_multi_chain_fields
    first = multi_state.states[0]
    levels = first.tempering_levels
    C = multi_state.num_chains
    c = jnp.arange(C)
    level = c % levels
    E = multi_chain_energy(multi_state)
    beta = jnp.asarray(first.temperature, dtype=jnp.float32)
    parity = sweep_index % 2
    is_lead = (level % 2 == parity) & (level < levels - 1)
    nxt = jnp.minimum(c + 1, C - 1)
    log_ratio = (beta - beta[nxt]) * (E[nxt] - E)
    u = random.uniform(key, (C,), dtype=jnp.float32)
    accept_lead = is_lead & (jnp.log(u) < log_ratio)
    accept_follow = jnp.concatenate([jnp.zeros(1, bool), accept_lead[:-1]])
    perm = jnp.where(accept_lead, nxt, c)
    perm = jnp.where(accept_follow, c - 1, perm)
    per, shared = split_multi_chain_fields(multi_state)
    per = jax.tree_util.tree_map(lambda a: a[perm], per)
    new_state = eqx.combine(per, shared)
    return multi_set_temperature(new_state, beta), accept_lead
