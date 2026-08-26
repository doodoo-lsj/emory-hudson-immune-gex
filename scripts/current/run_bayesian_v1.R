suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source("R/pca_mediation_pipeline.R")
source("R/bayesian_mediation_v1.R")

run_bayesian_v1 <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    coord_x = "imagecol",
    coord_y = "imagerow",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    standardize_covariates = TRUE,
    standardize_mediators = FALSE,
    adapt_delta = NULL,
    max_treedepth = NULL,
    use_cmdstanr_if_available = TRUE,
    posterior_mc_draws = NULL,
    refresh = 100) {
  total_start_time <- Sys.time()

  total_elapsed <- system.time({
    check_bayesian_v1_packages(require_brms = TRUE)

    backend <- NULL
    if (isTRUE(use_cmdstanr_if_available)) {
      if (requireNamespace("cmdstanr", quietly = TRUE)) {
        backend <- "cmdstanr"
        message("Using brms backend = 'cmdstanr'.")
      } else {
        message("cmdstanr is not installed; using brms default backend.")
      }
    }

    analysis_input <- prepare_bayesian_mediation_v1_data(
      rds_path = rds_path,
      coord_x = coord_x,
      coord_y = coord_y,
      standardize_covariates = standardize_covariates,
      standardize_mediators = standardize_mediators
    )

    dat_full_allpc <- analysis_input$dat_full_allpc
    dat_bayes <- analysis_input$dat_bayes
    formulas <- bayesian_mediation_v1_formulas(
      standardize_covariates = standardize_covariates
    )

    cat("\nDeterministic PCA validation\n")
    print(analysis_input$pca_validation)

    cat("\nPredictor scales in current complete-case analysis data\n")
    print(analysis_input$predictor_scale_summary)

    if (isTRUE(standardize_covariates)) {
      cat("\nCovariate standardization used for Bayesian fitting\n")
      print(analysis_input$covariate_scale_info)
      cat("\nStandardized covariate validation\n")
      print(analysis_input$standardized_covariate_validation)
    } else {
      cat("\nCovariate standardization disabled; Bayesian fitting uses raw x_coord/y_coord.\n")
    }

    priors <- make_bayesian_mediation_v1_priors(
      data = dat_bayes,
      standardize_covariates = standardize_covariates
    )

    frequentist <- compute_bayes_v1_frequentist_estimates(dat_full_allpc)

    cat("\nCurrent frequentist coordinate-adjusted decomposition\n")
    print(
      frequentist$summary |>
        select(contrast, TE, direct_effect = NDE, NIE, PM)
    )

    cat("\nCurrent frequentist mediator-specific paths\n")
    print(
      frequentist$path |>
        select(contrast, mediator, alpha = alpha_X_to_M, beta = beta_M_to_Y, indirect_component)
    )

    fits <- fit_bayesian_mediation_v1(
      data = dat_bayes,
      formulas = formulas,
      priors = priors,
      chains = chains,
      iter = iter,
      warmup = warmup,
      cores = cores,
      seed = seed,
      backend = backend,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      refresh = refresh
    )
    fit_timing <- attr(fits, "timing")

    posterior_processing_elapsed <- system.time({
      posterior_quantities <- compute_bayesian_mediation_v1_quantities(
        fits = fits,
        n_draws = posterior_mc_draws,
        seed = seed + 100
      )

      mediator_residual_dependence <- summarize_mediator_residual_dependence_v1(
        fits$mediator_joint
      )
    })

    path_summary <- posterior_quantities$path_summary
    decomposition_summary <- posterior_quantities$decomposition_summary

    frequentist_comparison <- compare_frequentist_bayesian_v1(
      frequentist = frequentist,
      path_summary = path_summary,
      decomposition_summary = decomposition_summary
    )

    diagnostics <- summarize_bayesian_v1_diagnostics(
      fits = fits,
      max_treedepth = max_treedepth %||% 10
    )
  })

  total_end_time <- Sys.time()

  timing <- dplyr::bind_rows(
    fit_timing,
    tibble::tibble(
      component = "posterior_draw_processing",
      elapsed_seconds = unname(posterior_processing_elapsed[["elapsed"]])
    ),
    tibble::tibble(
      component = "total_bayesian_analysis",
      elapsed_seconds = unname(total_elapsed[["elapsed"]])
    )
  )

  bayes_v1_results <- list(
    data = list(
      dat_full_allpc = dat_full_allpc,
      dat_bayes = dat_bayes,
      pca_df = analysis_input$pca_df,
      pca_fit = analysis_input$pca$pca_fit,
      predictor_scale_summary = analysis_input$predictor_scale_summary,
      covariate_scale_info = analysis_input$covariate_scale_info,
      pca_validation = analysis_input$pca_validation,
      standardized_covariate_validation = analysis_input$standardized_covariate_validation
    ),
    settings = list(
      rds_path = rds_path,
      coord_x = coord_x,
      coord_y = coord_y,
      chains = chains,
      iter = iter,
      warmup = warmup,
      cores = cores,
      seed = seed,
      standardize_covariates = standardize_covariates,
      standardize_mediators = standardize_mediators,
      backend = backend,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      posterior_mc_draws = posterior_mc_draws,
      formulas = formulas,
      priors = priors,
      start_time = total_start_time,
      end_time = total_end_time
    ),
    models = list(
      mediator_joint = fits$mediator_joint,
      outcome = fits$outcome,
      total = fits$total
    ),
    mediator_residual_dependence = mediator_residual_dependence,
    path_summary = path_summary,
    decomposition_summary = decomposition_summary,
    posterior_draws = list(
      tidy = posterior_quantities$draws_tidy,
      wide = posterior_quantities$draws_wide,
      path = posterior_quantities$path_draws,
      decomposition = posterior_quantities$decomposition_draws,
      sampled_indices = posterior_quantities$sampled_indices,
      available_draws = posterior_quantities$available_draws
    ),
    frequentist = frequentist,
    frequentist_comparison = frequentist_comparison,
    diagnostics = diagnostics,
    pm_diagnostics = posterior_quantities$pm_diagnostics,
    timing = timing
  )

  cat("\nMediator residual correlation summary\n")
  print(mediator_residual_dependence$correlation_summary)

  cat("\nPosterior mean residual correlation matrix\n")
  print(mediator_residual_dependence$correlation_mean_matrix)

  cat("\nPosterior mean residual covariance matrix\n")
  print(mediator_residual_dependence$covariance_mean_matrix)

  cat("\nBayesian v1 decomposition summary\n")
  print(decomposition_summary)

  cat("\nBayesian v1 mediator-specific path summary\n")
  print(path_summary)

  cat("\nPM diagnostics\n")
  print(posterior_quantities$pm_diagnostics)

  cat("\nFrequentist vs Bayesian v1 decomposition comparison\n")
  print(frequentist_comparison$decomposition)

  cat("\nFrequentist vs Bayesian v1 path comparison\n")
  print(frequentist_comparison$path)

  cat("\nMCMC diagnostics\n")
  print(diagnostics)

  cat("\nTiming\n")
  print(timing)

  cat("\nPrimary result object: bayes_v1_results\n")
  cat("Posterior predictive check examples to run manually:\n")
  cat("  brms::pp_check(bayes_v1_results$models$mediator_joint, resp = \"PC1R\")\n")
  cat("  brms::pp_check(bayes_v1_results$models$outcome)\n")
  cat("  brms::pp_check(bayes_v1_results$models$total)\n")

  bayes_v1_results
}

## Rscript scripts/current/run_bayesian_v1.R executes sampling.
## source("scripts/current/run_bayesian_v1.R") only defines run_bayesian_v1().
if (sys.nframe() == 0) {
  bayes_v1_results <- run_bayesian_v1()
}
