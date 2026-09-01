suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
})

source("R/pca_mediation_pipeline.R")
source("R/bayesian_mediation_v1.R")
source("R/sensitivity_v2.R")
source("R/sensitivity_A_residual_rho.R")
source("R/sensitivity_B_latent_U.R")
source("R/bayesian_sensitivity_A_B12.R")

run_bayesian_sensitivity <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    coord_x = "imagecol",
    coord_y = "imagerow",
    standardize_covariates = TRUE,
    chains = 4,
    iter = 1000,
    warmup = 500,
    cores = 4,
    seed = 123,
    adapt_delta = 0.9,
    max_treedepth = 10,
    refresh = 100,
    posterior_mc_draws = NULL,
    A_rho_grid = c(-0.3, 0, 0.3),
    B_p_U = 0.3,
    B_OR_UX_grid = c(1, 1.5),
    B_s_Y_grid = c(0, 0.3),
    B_s_M_grid = c(0, 0.3),
    run_B_sampling = FALSE,
    force_recompile_stan = FALSE) {
  check_bayesian_sensitivity_packages(require_brms = TRUE, require_stan = run_B_sampling)

  analysis_input <- prepare_bayesian_mediation_v1_data(
    rds_path = rds_path,
    coord_x = coord_x,
    coord_y = coord_y,
    standardize_covariates = standardize_covariates,
    standardize_mediators = FALSE
  )
  dat_bayes <- analysis_input$dat_bayes
  dat_full_allpc <- analysis_input$dat_full_allpc

  formulas <- bayesian_mediation_v1_formulas(standardize_covariates = standardize_covariates)
  priors <- make_bayesian_mediation_v1_priors(dat_bayes, standardize_covariates = standardize_covariates)

  cat("\nFitting baseline Bayesian mediation v1 models\n")
  baseline_fits <- fit_bayesian_mediation_v1(
    data = dat_bayes,
    formulas = formulas,
    priors = priors,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    backend = NULL,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = refresh
  )

  A_suite <- compute_bayesian_sensitivity_A_suite(
    fits = baseline_fits,
    n_draws = posterior_mc_draws,
    seed = seed + 100,
    rho_grid = A_rho_grid
  )
  baseline_quantities <- A_suite$baseline_quantities
  baseline_diagnostics <- summarize_bayesian_v1_diagnostics(
    fits = baseline_fits,
    max_treedepth = max_treedepth
  )

  cat("\nBayesian baseline decomposition, primary TE = reduced-form total-effect model coefficient\n")
  print(baseline_quantities$decomposition_summary)

  cat("\nBayesian baseline path summary\n")
  print(baseline_quantities$path_summary)

  cat("\nRunning Bayesian Sensitivity A posterior rho calculation\n")
  A_common <- A_suite$common
  A_one <- A_suite$one_at_a_time
  A_validation <- A_suite$validation

  cat("\nBayesian Sensitivity A rho=0 draw-by-draw validation\n")
  print(tibble::tibble(
    check = c("common vs baseline", "one-at-a-time vs baseline", "common vs one-at-a-time"),
    passed = c(
      A_validation$common_rho0_vs_baseline$passed,
      A_validation$one_at_a_time_rho0_vs_baseline$passed,
      A_validation$all_A_rho0_scenarios_identical$passed
    ),
    max_abs_draw_diff = c(
      A_validation$common_rho0_vs_baseline$max_abs_draw_diff,
      A_validation$one_at_a_time_rho0_vs_baseline$max_abs_draw_diff,
      A_validation$all_A_rho0_scenarios_identical$max_abs_draw_diff
    )
  ))

  cat("\nBayesian Sensitivity A common-rho summaries\n")
  print(A_common$summary)

  cat("\nBayesian Sensitivity A one-at-a-time summaries\n")
  print(A_one$summary)

  B1_results <- NULL
  B2_results <- NULL
  comparison <- list()

  if (isTRUE(run_B_sampling)) {
    stan_model <- bayes_sens_compile_B_model(force_recompile = force_recompile_stan)

    B1_grid <- expand.grid(
      p_U = B_p_U,
      OR_UX = B_OR_UX_grid,
      s_Y = B_s_Y_grid
    ) |>
      tibble::as_tibble()
    B2_grid <- expand.grid(
      p_U = B_p_U,
      OR_UX = B_OR_UX_grid,
      s_M = B_s_M_grid
    ) |>
      tibble::as_tibble()

    cat("\nRunning Bayesian Sensitivity B1 small grid\n")
    B1_raw <- purrr::pmap(B1_grid, function(p_U, OR_UX, s_Y) {
      fit_b <- fit_bayesian_B_scenario(
        data = dat_bayes,
        scenario = "B1_XY",
        p_U = p_U,
        OR_U_adj = OR_UX,
        OR_U_far = OR_UX,
        s_Y = s_Y,
        standardize_covariates = standardize_covariates,
        chains = chains,
        iter = iter,
        warmup = warmup,
        cores = cores,
        seed = seed + 400 + as.integer(round(100 * OR_UX)) + as.integer(round(100 * s_Y)),
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        refresh = refresh,
        stan_model = stan_model
      )
      quantities <- compute_bayesian_B_quantities(fit_b, n_draws = posterior_mc_draws, seed = seed + 500)
      list(fit = fit_b, quantities = quantities, diagnostics = summarize_rstan_sensitivity_diagnostics(fit_b$fit, max_treedepth = max_treedepth))
    })

    cat("\nRunning Bayesian Sensitivity B2 small grid\n")
    B2_raw <- purrr::pmap(B2_grid, function(p_U, OR_UX, s_M) {
      fit_b <- fit_bayesian_B_scenario(
        data = dat_bayes,
        scenario = "B2_XM",
        p_U = p_U,
        OR_U_adj = OR_UX,
        OR_U_far = OR_UX,
        s_M = s_M,
        standardize_covariates = standardize_covariates,
        chains = chains,
        iter = iter,
        warmup = warmup,
        cores = cores,
        seed = seed + 600 + as.integer(round(100 * OR_UX)) + as.integer(round(100 * s_M)),
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        refresh = refresh,
        stan_model = stan_model
      )
      quantities <- compute_bayesian_B_quantities(fit_b, n_draws = posterior_mc_draws, seed = seed + 700)
      list(fit = fit_b, quantities = quantities, diagnostics = summarize_rstan_sensitivity_diagnostics(fit_b$fit, max_treedepth = max_treedepth))
    })

    B1_results <- list(
      grid = B1_grid,
      raw = B1_raw,
      summary = bind_rows(map(B1_raw, ~ .x$quantities$summary)),
      path = bind_rows(map(B1_raw, ~ .x$quantities$path)),
      diagnostics = bind_rows(map(B1_raw, "diagnostics"), .id = "grid_id"),
      validation = do.call(
        validate_bayes_sens_TE_identity,
        setNames(map(B1_raw, "quantities"), paste0("B1_", seq_along(B1_raw)))
      )
    )
    B2_results <- list(
      grid = B2_grid,
      raw = B2_raw,
      summary = bind_rows(map(B2_raw, ~ .x$quantities$summary)),
      path = bind_rows(map(B2_raw, ~ .x$quantities$path)),
      diagnostics = bind_rows(map(B2_raw, "diagnostics"), .id = "grid_id"),
      validation = do.call(
        validate_bayes_sens_TE_identity,
        setNames(map(B2_raw, "quantities"), paste0("B2_", seq_along(B2_raw)))
      )
    )

    cat("\nB1 diagnostics\n")
    print(B1_results$diagnostics)
    cat("\nB2 diagnostics\n")
    print(B2_results$diagnostics)
  } else {
    message("run_B_sampling = FALSE: B1/B2 Stan models were not sampled. Set TRUE for the small diagnostic B grid.")
  }

  frequentist_A <- tryCatch(
    run_sensitivity_A_residual_rho(
      data = dat_full_allpc,
      rho_grid = A_rho_grid
    ),
    error = function(e) e
  )

  bayesian_sensitivity_results <- list(
    baseline = list(
      fits = baseline_fits,
      quantities = baseline_quantities,
      diagnostics = baseline_diagnostics
    ),
    A = list(
      common = A_common,
      one_at_a_time = A_one,
      validation = A_validation
    ),
    B1_XY = B1_results,
    B2_XM = B2_results,
    comparison = list(
      frequentist_A = frequentist_A,
      note = "Frequentist-vs-Bayesian comparison is qualitative; exact equality is not used as a validation criterion."
    ),
    data = list(
      dat_bayes = dat_bayes,
      dat_full_allpc = dat_full_allpc,
      predictor_scale_summary = analysis_input$predictor_scale_summary,
      covariate_scale_info = analysis_input$covariate_scale_info
    ),
    settings = list(
      rds_path = rds_path,
      coord_x = coord_x,
      coord_y = coord_y,
      standardize_covariates = standardize_covariates,
      chains = chains,
      iter = iter,
      warmup = warmup,
      cores = cores,
      seed = seed,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      posterior_mc_draws = posterior_mc_draws,
      A_rho_grid = A_rho_grid,
      B_p_U = B_p_U,
      B_OR_UX_grid = B_OR_UX_grid,
      B_s_Y_grid = B_s_Y_grid,
      B_s_M_grid = B_s_M_grid,
      run_B_sampling = run_B_sampling
    )
  )

  cat("\nA rho=0 feasibility/common summaries available in bayesian_sensitivity_results$A\n")
  cat("\nPrimary result object: bayesian_sensitivity_results\n")
  bayesian_sensitivity_results
}

## Small diagnostic example:
## bayesian_sensitivity_results <- run_bayesian_sensitivity(
##   chains = 2,
##   iter = 500,
##   warmup = 250,
##   cores = 2,
##   posterior_mc_draws = 500,
##   A_rho_grid = c(-0.3, 0, 0.3),
##   B_OR_UX_grid = c(1, 1.5),
##   B_s_Y_grid = c(0, 0.3),
##   B_s_M_grid = c(0, 0.3),
##   run_B_sampling = TRUE
## )

if (sys.nframe() == 0) {
  bayesian_sensitivity_results <- run_bayesian_sensitivity()
}
