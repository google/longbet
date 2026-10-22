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

# The Python fit treats a RankKLongBetConfig field that differs from its default
# as an explicit choice (and warns when such a choice is ignored), so the R
# defaults must be the Python ones: NULL where the Python default is None (the
# variant-dependent fields, which RankKLongBetConfig resolves to the budget of
# the selected variant) and the same value everywhere else.

# Python RankKLongBetConfig fields that the R front doors do not expose.
rank_k_python_only_fields <- c(
  "baseline_cells", "between_var_threshold", "catt_floor", "lengthscale_grid_km",
  "num_gfr", "ridge_eps", "tau_clip", "tau_prior_mean", "tau_prior_sd",
  "theta_f_clip", "theta_prior_mean", "theta_prior_sd", "warm_burnin",
  # Remediation options (ATT prior families, short-cohort policy, sampler budget and
  # treatment-forest prior controls, exposure-basis generators, projection accounting)
  "tau_prior_family", "short_cohort_policy", "thin", "sigma_nu_total",
  "sigma_nu_dynamic_total", "alpha_nu", "beta_nu", "nu_working_zeta_weight",
  "learn_sigma_nu", "sigma_nu_prior_scale", "treated_extra_variance", "basis_spec",
  "keep_unprojected_draws"
)

# Named list: field name -> default of RankKLongBetConfig (NULL for None). The
# dataclasses module is imported with conversion, so the tuple of fields arrives
# as an R list and each default as an R value.
rank_k_python_defaults <- function() {
  fields <- reticulate::import("dataclasses")$fields(
    reticulate::import("longbet")$RankKLongBetConfig
  )
  stats::setNames(lapply(fields, function(f) f$default),
                  vapply(fields, function(f) f$name, character(1)))
}

# A resolved RankKLongBetConfig as a named R list (NULL for None), as in
# R/encourage_model.R.
rank_k_config_list <- function(cfg) {
  reticulate::import("dataclasses")$asdict(cfg)
}

rank_k_budget_names <- c("engine", "uncertainty_propagation", "num_chains", "num_burnin",
                         "num_trees_mu", "num_trees_tau_mode1", "num_trees_tau_higher")

test_that("Rank-K R defaults mirror RankKLongBetConfig()", {
  skip_without_engine()
  py_defaults <- rank_k_python_defaults()
  config_formals <- as.list(formals(longbet:::.rank_k_config))
  fit_formals <- as.list(formals(longbet::fit_rank_k_longbet))
  config_names <- names(config_formals)

  # Every R configuration argument is a Python field, and only the documented
  # fields are Python-only.
  expect_true(all(config_names %in% names(py_defaults)))
  expect_setequal(setdiff(names(py_defaults), config_names), rank_k_python_only_fields)
  # fit_rank_k_longbet() exposes every configuration argument with the same default.
  expect_true(all(config_names %in% names(fit_formals)))
  expect_identical(fit_formals[config_names], config_formals)

  for (nm in config_names) {
    r_default <- eval(config_formals[[nm]])
    py_default <- py_defaults[[nm]]
    info <- sprintf("default of `%s`", nm)
    if (is.null(py_default)) {
      expect_null(r_default, info = info)
    } else {
      expect_equal(r_default, py_default, info = info)
    }
  }
  expect_true(all(c("engine", "uncertainty_propagation", "max_forecast_draws", "num_chains",
                    "num_burnin") %in% config_names))
})

test_that(".rank_k_config() resolves each variant like RankKLongBetConfig()", {
  skip_without_engine()
  lb <- reticulate::import("longbet", convert = FALSE)

  propagated <- rank_k_config_list(longbet:::.rank_k_config())
  expect_equal(propagated, rank_k_config_list(lb$RankKLongBetConfig()))
  expect_equal(propagated[rank_k_budget_names],
               list(engine = "jax", uncertainty_propagation = TRUE, num_chains = 2L,
                    num_burnin = 100L, num_trees_mu = 50L, num_trees_tau_mode1 = 30L,
                    num_trees_tau_higher = 15L))
  expect_equal(propagated[c("n_forecast_draws", "draws_per_component")],
               list(n_forecast_draws = 8L, draws_per_component = 25L))
  expect_null(propagated$max_forecast_draws)

  plug_in <- rank_k_config_list(longbet:::.rank_k_config(uncertainty_propagation = FALSE))
  expect_equal(plug_in, rank_k_config_list(lb$RankKLongBetConfig(uncertainty_propagation = FALSE)))
  expect_equal(plug_in[rank_k_budget_names],
               list(engine = "stochtree", uncertainty_propagation = FALSE, num_chains = 1L,
                    num_burnin = 200L, num_trees_mu = 75L, num_trees_tau_mode1 = 35L,
                    num_trees_tau_higher = 20L))
  expect_equal(plug_in$num_mcmc, 600L)
  expect_null(plug_in$max_forecast_draws)

  # `engine` alone still selects its variant (backward compatibility).
  expect_equal(rank_k_config_list(longbet:::.rank_k_config(engine = "jax")), propagated)
  expect_equal(rank_k_config_list(longbet:::.rank_k_config(uncertainty_propagation = TRUE)),
               propagated)
  expect_equal(rank_k_config_list(longbet:::.rank_k_config(engine = "stochtree")), plug_in)
})

test_that("Explicit Rank-K budgets override the variant defaults", {
  skip_without_engine()
  application <- rank_k_config_list(longbet:::.rank_k_config(
    num_chains = 4L, num_burnin = 300L, draws_per_component = 100L
  ))
  expect_equal(application[c("engine", "num_chains", "num_burnin", "draws_per_component",
                             "num_trees_mu")],
               list(engine = "jax", num_chains = 4L, num_burnin = 300L,
                    draws_per_component = 100L, num_trees_mu = 50L))
  expanded <- rank_k_config_list(longbet:::.rank_k_config(max_forecast_draws = 16L))
  expect_equal(expanded$max_forecast_draws, 16L)
  plug_in <- rank_k_config_list(longbet:::.rank_k_config(uncertainty_propagation = FALSE,
                                                         num_burnin = 400L, num_trees_mu = 100L))
  expect_equal(plug_in[c("engine", "num_burnin", "num_trees_mu", "num_trees_tau_mode1")],
               list(engine = "stochtree", num_burnin = 400L, num_trees_mu = 100L,
                    num_trees_tau_mode1 = 35L))
})

test_that("Contradictory Rank-K variants and invalid flags are errors", {
  skip_without_engine()
  expect_error(longbet:::.rank_k_config(engine = "jax", uncertainty_propagation = FALSE),
               "not a supported variant")
  expect_error(longbet:::.rank_k_config(engine = "stochtree", uncertainty_propagation = TRUE),
               "not implemented")
  expect_error(longbet:::.rank_k_config(engine = "bartz"), "engine must be None")
  expect_error(longbet:::.rank_k_config(n_forecast_draws = 4L, max_forecast_draws = 2L),
               "max_forecast_draws")
  expect_error(longbet:::.rank_k_config(uncertainty_propagation = NA), "must be NULL, TRUE or FALSE")
  expect_error(longbet:::.rank_k_config(uncertainty_propagation = 1), "must be NULL, TRUE or FALSE")
  expect_error(longbet:::.rank_k_config(use_forecast = NA), "must be TRUE or FALSE")
  expect_error(longbet:::.rank_k_config(use_spatial_gp = c(TRUE, FALSE)), "must be TRUE or FALSE")

  # The front doors build the configuration before any computation.
  y <- matrix(0, 4, 4)
  x <- matrix(0, 4, 1)
  expect_error(fit_rank_k_longbet(y, y, x, engine = "jax", uncertainty_propagation = FALSE),
               "not a supported variant")
  expect_error(rank_k_placebo_audit(y, y, x, engine = "stochtree", uncertainty_propagation = TRUE),
               "not implemented")
})

test_that("Rank-K counts must be single whole numbers", {
  skip_without_engine()
  # Checked in R, because as.integer() would truncate 2.5 or pass NA on.
  expect_error(longbet:::.rank_k_config(rank_k = 2.5), "`rank_k` must be a single whole number")
  expect_error(longbet:::.rank_k_config(num_mcmc = "600"), "`num_mcmc` must be a single whole number")
  expect_error(longbet:::.rank_k_config(random_seed = c(1, 2)), "`random_seed` must be a single whole number")
  expect_error(longbet:::.rank_k_config(num_chains = NA), "`num_chains` must be NULL or a single whole number")
  expect_error(longbet:::.rank_k_config(num_trees_mu = TRUE),
               "`num_trees_mu` must be NULL or a single whole number")
  expect_error(longbet:::.rank_k_config(num_burnin = Inf), "`num_burnin` must be NULL or a single whole number")
  # The ranges are checked by RankKLongBetConfig.
  expect_error(longbet:::.rank_k_config(rank_k = 0L), "rank_k must be >= 1")
  expect_error(longbet:::.rank_k_config(num_threads = 0L), "num_threads must be >= 1, or -1")
  # Whole doubles are accepted and passed on as integers; 0 forecast draws means none.
  cfg <- rank_k_config_list(longbet:::.rank_k_config(rank_k = 3, num_burnin = 150, num_threads = -1,
                                                     n_forecast_draws = 0))
  expect_equal(cfg[c("rank_k", "num_burnin", "num_threads", "n_forecast_draws")],
               list(rank_k = 3L, num_burnin = 150L, num_threads = -1L, n_forecast_draws = 0L))
  # The placebo window length is checked in the same way, before any computation.
  y <- matrix(0, 4, 4)
  expect_error(rank_k_placebo_audit(y, y, matrix(0, 4, 1), pseudo_anchor_offset = 2.5),
               "`pseudo_anchor_offset` must be a single whole number")
})

test_that("fit_rank_k_longbet() keeps the 0.4.0 argument positions", {
  v040 <- c("y", "z", "x", "features", "unit_features", "ps", "baseline_forecast",
            "baseline_forecast_draws", "coords", "cluster_ids", "t_0_panel", "rank_k", "engine",
            "use_forecast", "uncertainty_propagation", "n_forecast_draws", "use_spatial_gp",
            "use_factor_loadings", "n_factor_loadings", "num_chains", "num_burnin", "num_mcmc",
            "num_trees_mu", "num_trees_tau_mode1", "num_trees_tau_higher", "random_seed",
            "num_threads")
  fit_names <- names(formals(longbet::fit_rank_k_longbet))
  expect_identical(fit_names[seq_along(v040)], v040)
  expect_setequal(fit_names[-seq_along(v040)],
                  c("pre_anchor_loadings", "cohorts", "max_forecast_draws", "draws_per_component"))
})
