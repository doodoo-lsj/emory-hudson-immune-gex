suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source("R/pca_mediation_pipeline.R")
source("R/bayesian_mediation_v0.R")

run_bayesian_v0 <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    coord_x = "imagecol",
    coord_y = "imagerow",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    standardize_coordinates = TRUE,
    adapt_delta = NULL,
    max_treedepth = NULL,
    use_cmdstanr_if_available = TRUE,
    posterior_mc_draws = NULL,
    refresh = 100) {
  total_start_time <- Sys.time()
  total_elapsed <- system.time({
    check_bayesian_v0_packages()

    backend <- NULL
    if (isTRUE(use_cmdstanr_if_available)) {
      if (requireNamespace("cmdstanr", quietly = TRUE)) {
        backend <- "cmdstanr"
        message("Using brms backend = 'cmdstanr'.")
      } else {
        message("cmdstanr is not installed; using brms default backend.")
      }
    }

    analysis_input <- prepare_bayesian_mediation_v0_data(
      rds_path = rds_path,
      coord_x = coord_x,
      coord_y = coord_y,
      standardize_coordinates = standardize_coordinates
    )

    dat_full_allpc <- analysis_input$dat_full_allpc
    dat_bayes <- analysis_input$dat_bayes
    formulas <- bayesian_mediation_v0_formulas(
      use_standardized_coordinates = standardize_coordinates
    )

    cat("\nPredictor scales in the current analysis data\n")
    print(analysis_input$predictor_scale_summary)

    if (isTRUE(standardize_coordinates)) {
      cat("\nCoordinate standardization used for Bayesian fitting\n")
      print(analysis_input$coordinate_scale_info)
      cat(
        "\nBayesian fitting uses x_coord_std/y_coord_std for numerical geometry. ",
        "Exposure coefficients and PC mediator coefficients remain on the original mediation scale.\n",
        sep = ""
      )
    } else {
      cat("\nCoordinate standardization is disabled; Bayesian fitting uses original x_coord/y_coord.\n")
    }

    # Priors are generated on the actual fitting scale. With standardized
    # coordinates, the coordinate priors are placed on x_coord_std/y_coord_std,
    # which improves sampling geometry without changing the regression
    # estimand. To modify priors for v0, edit this object before fitting.
    priors <- make_bayesian_mediation_v0_priors(
      data = dat_bayes,
      formulas = formulas
    )

    frequentist <- compute_frequentist_current_estimates(dat_full_allpc)

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

    # Divergences were not observed in the previous v0 run; the primary change
    # here is coordinate centering/scaling. Keep sampler controls at brms
    # defaults unless diagnostics after rerun still show treedepth/Rhat issues.
    fits <- fit_bayesian_mediation_v0(
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
      posterior_quantities <- compute_bayesian_mediation_quantities(
        fits = fits,
        n_draws = posterior_mc_draws,
        seed = seed + 100
      )
      fixed_effect_summary_original_scale <- summarize_original_scale_fixed_effects(
        fits = fits,
        coordinate_scale_info = analysis_input$coordinate_scale_info
      )
    })

    path_summary <- posterior_quantities$path_summary
    decomposition_summary <- posterior_quantities$decomposition_summary

    frequentist_comparison <- compare_frequentist_bayesian_v0(
      frequentist = frequentist,
      bayesian_path_summary = path_summary,
      bayesian_decomposition_summary = decomposition_summary
    )

    diagnostics <- summarize_bayesian_diagnostics(
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

  bayes_v0_results <- list(
    config = list(
      rds_path = rds_path,
      coord_x = coord_x,
      coord_y = coord_y,
      chains = chains,
      iter = iter,
      warmup = warmup,
      cores = cores,
      seed = seed,
      standardize_coordinates = standardize_coordinates,
      backend = backend,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      posterior_mc_draws = posterior_mc_draws,
      formulas = formulas,
      priors = priors,
      start_time = total_start_time,
      end_time = total_end_time
    ),
    data = list(
      dat_full_allpc = dat_full_allpc,
      dat_bayes = dat_bayes,
      pca_df = analysis_input$pca_df,
      pca_fit = analysis_input$pca$pca_fit,
      predictor_scale_summary = analysis_input$predictor_scale_summary,
      coordinate_scale_info = analysis_input$coordinate_scale_info
    ),
    frequentist = frequentist,
    fits = fits,
    posterior_draws = list(
      mediation_quantities = posterior_quantities$draws,
      sampled_indices = posterior_quantities$sampled_indices,
      available_draws = posterior_quantities$available_draws
    ),
    path_summary = path_summary,
    decomposition_summary = decomposition_summary,
    fixed_effect_summary_original_scale = fixed_effect_summary_original_scale,
    frequentist_comparison = frequentist_comparison,
    diagnostics = diagnostics,
    timing = timing
  )

  cat("\nBayesian decomposition summary\n")
  print(decomposition_summary)

  cat("\nBayesian path summary\n")
  print(path_summary)

  cat("\nFrequentist vs Bayesian decomposition comparison\n")
  print(frequentist_comparison$decomposition)

  cat("\nFrequentist vs Bayesian path comparison\n")
  print(frequentist_comparison$path)

  cat("\nMCMC diagnostics\n")
  print(diagnostics)

  cat("\nTiming\n")
  print(timing)

  cat("\nPrimary result object: bayes_v0_results\n")
  cat("Posterior predictive check examples to run manually:\n")
  cat("  brms::pp_check(bayes_v0_results$fits$outcome)\n")
  cat("  brms::pp_check(bayes_v0_results$fits$total)\n")
  cat("  brms::pp_check(bayes_v0_results$fits$mediator_PC1_R)\n")

  bayes_v0_results
}

## Rscript scripts/current/run_bayesian_v0.R executes sampling.
## source("scripts/current/run_bayesian_v0.R") only defines run_bayesian_v0().
if (sys.nframe() == 0) {
  bayes_v0_results <- run_bayesian_v0()
}
