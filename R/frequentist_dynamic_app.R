# Dynamic frequentist mediation runtime helpers for prepared Shiny data.
#
# This layer mirrors the Shiny Bayesian runtime data specification, but uses the
# existing lm-based frequentist decomposition and residual-rho sensitivity.

frequentist_dynamic_model_version <- function() {
  "frequentist_mediation_dynamic_v1"
}

frequentist_dynamic_prepare_data <- function(analysis_obj, selected_pcs, pc_signs) {
  prepared <- bayesian_dynamic_prepare_data(
    analysis_obj = analysis_obj,
    selected_pcs = selected_pcs,
    pc_signs = pc_signs,
    standardize_numeric_covariates = TRUE
  )
  prepared$model_version <- frequentist_dynamic_model_version()
  prepared
}

map_frequentist_mediator_names <- function(df, mediator_map, columns = c("mediator", "varied_mediator")) {
  name_map <- stats::setNames(mediator_map$mediator, mediator_map$model_name)
  out <- df
  for (column in intersect(columns, names(out))) {
    mapped <- unname(name_map[out[[column]]])
    out[[column]] <- ifelse(is.na(mapped), out[[column]], mapped)
  }
  out
}

make_frequentist_dynamic_fit <- function(analysis_obj, selected_pcs, pc_signs) {
  prepared <- frequentist_dynamic_prepare_data(analysis_obj, selected_pcs, pc_signs)
  fit <- decompose_linear_multix(
    data = prepared$data,
    exposures = prepared$exposures,
    outcome = "Y",
    mediators = prepared$mediator_model_names,
    covariates = prepared$covariate_terms
  )

  fit$summary <- fit$summary |>
    dplyr::mutate(contrast = unname(prepared$exposure_labels[exposure]))
  fit$summary$contrast <- ifelse(is.na(fit$summary$contrast), fit$summary$exposure, fit$summary$contrast)
  fit$path <- fit$path |>
    map_frequentist_mediator_names(prepared$mediator_map) |>
    dplyr::mutate(contrast = unname(prepared$exposure_labels[exposure]))
  fit$path$contrast <- ifelse(is.na(fit$path$contrast), fit$path$exposure, fit$path$contrast)

  list(
    source = "dynamic_fit",
    model_id = frequentist_dynamic_model_version(),
    prepared = prepared,
    fit = fit,
    selected_pcs = selected_pcs,
    created_at = Sys.time()
  )
}

make_frequentist_example_fit <- function(analysis_obj) {
  fit <- analysis_obj$mediation_coord
  fit$summary <- fit$summary |>
    dplyr::mutate(contrast = unname(analysis_obj$contrast_labels[exposure]))
  fit$summary$contrast <- ifelse(is.na(fit$summary$contrast), fit$summary$exposure, fit$summary$contrast)
  fit$path <- fit$path |>
    dplyr::mutate(
      mediator = dplyr::case_when(
        mediator == "PC1_R" ~ "PC1",
        mediator == "PC2_R" ~ "PC2",
        TRUE ~ mediator
      ),
      contrast = unname(analysis_obj$contrast_labels[exposure])
    )
  fit$path$contrast <- ifelse(is.na(fit$path$contrast), fit$path$exposure, fit$path$contrast)

  prepared <- frequentist_dynamic_prepare_data(
    analysis_obj,
    selected_pcs = c("PC1", "PC2", "PC3"),
    pc_signs = stats::setNames(rep(1, 6), paste0("PC", seq_len(6)))
  )

  list(
    source = "validated_example",
    model_id = "frequentist_example_current",
    prepared = prepared,
    fit = fit,
    selected_pcs = c("PC1", "PC2", "PC3"),
    created_at = Sys.time()
  )
}

frequentist_dynamic_cache_key <- function(analysis_obj, selected_pcs, pc_signs) {
  digest::digest(
    list(
      model_version = frequentist_dynamic_model_version(),
      data_source = analysis_obj$spec$data_source,
      metadata = analysis_obj$spec$analysis_data$metadata,
      feature_matrix = analysis_obj$spec$analysis_data$feature_matrix,
      exposures = analysis_obj$exposures,
      covariates = analysis_obj$covariates,
      selected_pcs = selected_pcs,
      pc_signs = unname(pc_signs[selected_pcs])
    ),
    algo = "xxhash64"
  )
}

bootstrap_frequentist_dynamic <- function(frequentist_result,
                                          B = 1000,
                                          seed = 123,
                                          progress_callback = NULL) {
  prepared <- frequentist_result$prepared
  boot <- bootstrap_decomposition(
    data = prepared$data,
    B = B,
    seed = seed,
    exposures = prepared$exposures,
    outcome = "Y",
    mediators = prepared$mediator_model_names,
    covariates = prepared$covariate_terms,
    progress_callback = progress_callback
  )
  ci <- summarize_bootstrap_ci(boot$boot_summary_df, boot$boot_path_df)
  ci$summary_ci <- ci$summary_ci |>
    dplyr::mutate(contrast = unname(prepared$exposure_labels[exposure]))
  ci$summary_ci$contrast <- ifelse(is.na(ci$summary_ci$contrast), ci$summary_ci$exposure, ci$summary_ci$contrast)
  ci$path_ci <- ci$path_ci |>
    map_frequentist_mediator_names(prepared$mediator_map) |>
    dplyr::mutate(contrast = unname(prepared$exposure_labels[exposure]))
  ci$path_ci$contrast <- ifelse(is.na(ci$path_ci$contrast), ci$path_ci$exposure, ci$path_ci$contrast)

  list(bootstrap = boot, ci = ci)
}

compute_frequentist_dynamic_sensitivity <- function(frequentist_result,
                                                    rho_grid = seq(-0.3, 0.3, by = 0.05)) {
  prepared <- frequentist_result$prepared

  common <- residual_corr_sensitivity(
    data = prepared$data,
    exposures = prepared$exposures,
    outcome = "Y",
    mediators = prepared$mediator_model_names,
    covariates = prepared$covariate_terms,
    rho_grid = rho_grid,
    mode = "common"
  )
  one_at_a_time <- residual_corr_sensitivity(
    data = prepared$data,
    exposures = prepared$exposures,
    outcome = "Y",
    mediators = prepared$mediator_model_names,
    covariates = prepared$covariate_terms,
    rho_grid = rho_grid,
    mode = "one_at_a_time"
  )

  rho0_common <- validate_rho0_matches_observed(common)
  rho0_one_at_a_time <- validate_rho0_matches_observed(one_at_a_time)

  map_result <- function(x) {
    x$observed_path <- x$observed_path |>
      map_frequentist_mediator_names(prepared$mediator_map)
    x$sensitivity_summary <- x$sensitivity_summary |>
      map_frequentist_mediator_names(prepared$mediator_map, columns = "varied_mediator")
    x$sensitivity_path <- x$sensitivity_path |>
      map_frequentist_mediator_names(prepared$mediator_map)
    x$sensitivity_beta <- x$sensitivity_beta |>
      map_frequentist_mediator_names(prepared$mediator_map)
    x
  }

  common <- map_result(common)
  one_at_a_time <- map_result(one_at_a_time)

  list(
    common = common,
    one_at_a_time = one_at_a_time,
    rho0_common = rho0_common,
    rho0_one_at_a_time = rho0_one_at_a_time,
    tipping_common = summarize_rho_tipping_points(common),
    tipping_one_at_a_time = summarize_rho_tipping_points(one_at_a_time)
  )
}

validate_frequentist_sign_reversal <- function(analysis_obj, selected_pcs, pc_signs, pc_to_reverse) {
  if (!pc_to_reverse %in% selected_pcs) {
    stop("pc_to_reverse must be one of the selected PCs.", call. = FALSE)
  }
  base <- make_frequentist_dynamic_fit(analysis_obj, selected_pcs, pc_signs)
  flipped_signs <- pc_signs
  flipped_signs[[pc_to_reverse]] <- -flipped_signs[[pc_to_reverse]]
  flipped <- make_frequentist_dynamic_fit(analysis_obj, selected_pcs, flipped_signs)

  base_path <- base$fit$path
  flipped_path <- flipped$fit$path
  joined_path <- dplyr::left_join(
    base_path,
    flipped_path,
    by = c("exposure", "contrast", "mediator"),
    suffix = c("_base", "_flipped")
  )
  joined_summary <- dplyr::left_join(
    base$fit$summary,
    flipped$fit$summary,
    by = c("exposure", "contrast"),
    suffix = c("_base", "_flipped")
  )

  reversed_rows <- joined_path$mediator == pc_to_reverse
  list(
    alpha_flips = all(abs(joined_path$alpha_X_to_M_base[reversed_rows] + joined_path$alpha_X_to_M_flipped[reversed_rows]) < 1e-8),
    beta_flips = all(abs(joined_path$beta_M_to_Y_base[reversed_rows] + joined_path$beta_M_to_Y_flipped[reversed_rows]) < 1e-8),
    indirect_invariant = all(abs(joined_path$indirect_component_base - joined_path$indirect_component_flipped) < 1e-8),
    decomposition_invariant = all(abs(joined_summary$TE_base - joined_summary$TE_flipped) < 1e-8) &&
      all(abs(joined_summary$NDE_base - joined_summary$NDE_flipped) < 1e-8) &&
      all(abs(joined_summary$NIE_base - joined_summary$NIE_flipped) < 1e-8) &&
      all(abs(joined_summary$PM_base - joined_summary$PM_flipped) < 1e-8)
  )
}
