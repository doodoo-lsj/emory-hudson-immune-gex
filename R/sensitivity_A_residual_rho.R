# Sensitivity A: mediator-outcome residual-correlation sensitivity.
#
# This file preserves R/sensitivity_v2.R and wraps its current residual-rho
# logic into a clearer v3-style result object for the joint-mediator baseline.
# Mathematically, sensitivity_v2.R estimates the mediator residual covariance
# Sigma_M from the three current mediator regressions and perturbs beta by
# beta(rho) = beta_hat - Sigma_M^{-1} c_MY, where
# c_MY,j = rho_j * sd(error_Mj) * sd(error_Y).

check_sensitivity_A_packages <- function() {
  required <- c("dplyr", "tidyr", "tibble", "purrr", "ggplot2")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Sensitivity A requires installed R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

sensitivity_A_contrast_label <- function(exposure) {
  dplyr::case_when(
    exposure == "X_adjacent" ~ "adjacent vs inside",
    exposure == "X_far" ~ "far vs inside",
    TRUE ~ exposure
  )
}

prepare_current_sensitivity_data <- function(rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
                                             coord_x = "imagecol",
                                             coord_y = "imagerow") {
  needed <- c(
    "load_spatial_mediation_rds",
    "prepare_analysis_metadata",
    "fit_pca_mediators",
    "build_pca_score_data",
    "build_full_analysis_data"
  )
  missing <- needed[!vapply(needed, exists, logical(1), mode = "function")]
  if (length(missing) > 0) {
    stop("Missing current pipeline functions: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  gex <- load_spatial_mediation_rds(rds_path)
  analysis_df <- prepare_analysis_metadata(gex, coord_x = coord_x, coord_y = coord_y)
  pca <- fit_pca_mediators(gex$M_expr)
  pca_df <- build_pca_score_data(analysis_df, pca$pca_fit)
  dat_full_allpc <- build_full_analysis_data(pca_df)

  list(
    gex = gex,
    analysis_df = analysis_df,
    pca = pca,
    pca_df = pca_df,
    dat_full_allpc = dat_full_allpc
  )
}

fit_current_joint_baseline <- function(data,
                                       exposures = c("X_adjacent", "X_far"),
                                       outcome = "Y",
                                       mediators = c("PC1_R", "PC2_R", "PC3"),
                                       covariates = c("x_coord", "y_coord")) {
  if (!exists("decompose_linear_multix", mode = "function")) {
    stop("decompose_linear_multix() must be available from R/pca_mediation_pipeline.R.", call. = FALSE)
  }

  decomposition <- decompose_linear_multix(
    data = data,
    exposures = exposures,
    outcome = outcome,
    mediators = mediators,
    covariates = covariates
  )

  mediator_formula <- stats::as.formula(
    paste0("cbind(", paste(mediators, collapse = ", "), ") ~ ", paste(c(exposures, covariates), collapse = " + "))
  )
  joint_mediator_fit <- stats::lm(mediator_formula, data = data)
  mediator_residuals <- stats::resid(joint_mediator_fit)
  colnames(mediator_residuals) <- mediators

  list(
    decomposition = decomposition,
    joint_mediator_fit = joint_mediator_fit,
    Sigma_M_mle = crossprod(mediator_residuals) / nrow(mediator_residuals),
    Sigma_M_unbiased = stats::cov(mediator_residuals),
    residual_correlation = stats::cor(mediator_residuals),
    summary = decomposition$summary |>
      dplyr::mutate(contrast = sensitivity_A_contrast_label(exposure)),
    path = decomposition$path |>
      dplyr::mutate(contrast = sensitivity_A_contrast_label(exposure))
  )
}

add_sensitivity_A_changes <- function(result) {
  baseline_summary <- result$observed_summary |>
    dplyr::select(exposure, baseline_TE = TE, baseline_direct = NDE, baseline_NIE = NIE, baseline_PM = PM)
  baseline_path <- result$observed_path |>
    dplyr::select(
      exposure,
      mediator,
      baseline_alpha = alpha_X_to_M,
      baseline_beta = beta_M_to_Y,
      baseline_indirect = indirect_component
    )

  sensitivity_summary <- result$sensitivity_summary |>
    dplyr::mutate(contrast = sensitivity_A_contrast_label(exposure)) |>
    dplyr::left_join(baseline_summary, by = "exposure") |>
    dplyr::mutate(
      direct_effect = NDE_adjusted,
      total_IIE = NIE_adjusted,
      PM = PM_adjusted,
      total_IIE_change = total_IIE - baseline_NIE,
      total_IIE_relative_change = total_IIE_change / baseline_NIE,
      direct_change = direct_effect - baseline_direct,
      PM_change = PM - baseline_PM
    )

  sensitivity_path <- result$sensitivity_path |>
    dplyr::mutate(contrast = sensitivity_A_contrast_label(exposure)) |>
    dplyr::left_join(baseline_path, by = c("exposure", "mediator")) |>
    dplyr::mutate(
      alpha = alpha_X_to_M,
      beta = beta_adjusted,
      IIE = indirect_adjusted,
      beta_change = beta - baseline_beta,
      IIE_change = IIE - baseline_indirect,
      IIE_relative_change = IIE_change / baseline_indirect
    )

  list(summary = sensitivity_summary, path = sensitivity_path)
}

run_sensitivity_A_residual_rho <- function(data,
                                           rho_grid = seq(-0.5, 0.5, by = 0.05),
                                           exposures = c("X_adjacent", "X_far"),
                                           outcome = "Y",
                                           mediators = c("PC1_R", "PC2_R", "PC3"),
                                           covariates = c("x_coord", "y_coord"),
                                           tolerance = 1e-8) {
  check_sensitivity_A_packages()
  needed <- c("residual_corr_sensitivity", "validate_rho0_matches_observed", "summarize_rho_tipping_points")
  missing <- needed[!vapply(needed, exists, logical(1), mode = "function")]
  if (length(missing) > 0) {
    stop("Missing sensitivity_v2 functions: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  baseline <- fit_current_joint_baseline(
    data = data,
    exposures = exposures,
    outcome = outcome,
    mediators = mediators,
    covariates = covariates
  )

  common <- residual_corr_sensitivity(
    data = data,
    exposures = exposures,
    outcome = outcome,
    mediators = mediators,
    covariates = covariates,
    rho_grid = rho_grid,
    mode = "common"
  )
  one_at_a_time <- residual_corr_sensitivity(
    data = data,
    exposures = exposures,
    outcome = outcome,
    mediators = mediators,
    covariates = covariates,
    rho_grid = rho_grid,
    mode = "one_at_a_time"
  )

  common_changes <- add_sensitivity_A_changes(common)
  one_changes <- add_sensitivity_A_changes(one_at_a_time)

  validation <- list(
    common_rho0 = validate_rho0_matches_observed(common, tolerance = tolerance),
    one_at_a_time_rho0 = validate_rho0_matches_observed(one_at_a_time, tolerance = tolerance)
  )
  validation$passed <- isTRUE(validation$common_rho0$passed) && isTRUE(validation$one_at_a_time_rho0$passed)

  list(
    method_specification = list(
      current_sensitivity_v2_math = paste(
        "rho_j is interpreted as Corr(error_Mj, error_Y).",
        "Sigma_M is estimated as the full covariance matrix of mediator residuals.",
        "For fixed c_MY,j = rho_j * sd(error_Mj) * sd(error_Y), beta is perturbed as beta_hat - solve(Sigma_M, c_MY).",
        "Feasibility is checked using c_MY' Sigma_M^{-1} c_MY / sigma_Y^2 <= 1."
      ),
      equals_desired_A = TRUE,
      uncertainty_note = "This residual-rho grid is a deterministic sensitivity curve around fitted frequentist estimates; bootstrap CIs are not newly generated here."
    ),
    baseline = baseline,
    sensitivity_A_common_rho = list(
      raw = common,
      summary = common_changes$summary,
      path = common_changes$path,
      tipping = summarize_rho_tipping_points(common)
    ),
    sensitivity_A_one_at_a_time = list(
      raw = one_at_a_time,
      summary = one_changes$summary,
      path = one_changes$path,
      tipping = summarize_rho_tipping_points(one_at_a_time)
    ),
    validation = validation,
    plots = build_rho_sensitivity_plots(common, one_at_a_time)
  )
}
