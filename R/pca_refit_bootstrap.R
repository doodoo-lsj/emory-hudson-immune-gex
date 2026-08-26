# PCA-refit bootstrap for the current frequentist PCA mediation analysis.
#
# The existing fixed-PCA bootstrap in R/pca_mediation_pipeline.R is preserved.
# This helper adds a separate bootstrap that resamples spots, refits PCA on the
# bootstrap gene-expression matrix, aligns bootstrap PCs to the full-sample
# oriented reference PCs, and refits the current linear mediation decomposition.

check_pca_refit_bootstrap_packages <- function() {
  required <- c("dplyr", "tidyr", "tibble", "purrr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "PCA-refit bootstrap requires installed R packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

pca_refit_contrast_label <- function(exposure) {
  dplyr::case_when(
    exposure == "X_adjacent" ~ "adjacent vs inside",
    exposure == "X_far" ~ "far vs inside",
    TRUE ~ exposure
  )
}

prepare_pca_refit_bootstrap_reference <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    coord_x = "imagecol",
    coord_y = "imagerow",
    mediators = c("PC1_R", "PC2_R", "PC3")) {
  check_pca_refit_bootstrap_packages()

  needed_functions <- c(
    "load_spatial_mediation_rds",
    "prepare_analysis_metadata",
    "fit_pca_mediators",
    "build_pca_score_data",
    "build_loading_data",
    "build_full_analysis_data",
    "decompose_linear_multix",
    "bootstrap_decomposition"
  )
  missing_functions <- needed_functions[!vapply(needed_functions, exists, logical(1), mode = "function")]
  if (length(missing_functions) > 0) {
    stop(
      "Missing current pipeline functions: ",
      paste(missing_functions, collapse = ", "),
      ". Source R/pca_mediation_pipeline.R before running PCA-refit bootstrap.",
      call. = FALSE
    )
  }

  gex <- load_spatial_mediation_rds(rds_path)
  analysis_df <- prepare_analysis_metadata(gex, coord_x = coord_x, coord_y = coord_y)
  pca <- fit_pca_mediators(gex$M_expr)
  pca_df <- build_pca_score_data(analysis_df, pca$pca_fit)
  loading_df <- build_loading_data(pca$pca_fit)
  dat_full_allpc <- build_full_analysis_data(pca_df)

  pca_df_with_row_id <- pca_df |>
    dplyr::mutate(
      .row_id = dplyr::row_number(),
      X_cat = factor(X_cat, levels = c("inside", "adjacent", "far")),
      X_adjacent = as.integer(X_cat == "adjacent"),
      X_far = as.integer(X_cat == "far")
    ) |>
    dplyr::select(.row_id, Y, X_cat, X_adjacent, X_far, PC1_R, PC2_R, PC3, x_coord, y_coord) |>
    tidyr::drop_na()

  dat_with_row_id <- pca_df_with_row_id
  dat_without_row_id <- dat_with_row_id |>
    dplyr::select(-.row_id)

  if (!isTRUE(all.equal(dat_without_row_id, dat_full_allpc, check.attributes = FALSE))) {
    stop("Internal data reconstruction does not match build_full_analysis_data().", call. = FALSE)
  }

  M_expr_analysis <- gex$M_expr[dat_with_row_id$.row_id, , drop = FALSE]
  reference_loadings <- as.matrix(loading_df[, mediators, drop = FALSE])
  rownames(reference_loadings) <- loading_df$gene

  point_estimates <- decompose_linear_multix(
    data = dat_full_allpc,
    exposures = c("X_adjacent", "X_far"),
    outcome = "Y",
    mediators = mediators,
    covariates = c("x_coord", "y_coord")
  )

  list(
    gex = gex,
    analysis_df = analysis_df,
    pca = pca,
    pca_df = pca_df,
    loading_df = loading_df,
    dat_full_allpc = dat_full_allpc,
    dat_with_row_id = dat_with_row_id,
    M_expr_analysis = M_expr_analysis,
    reference_loadings = reference_loadings,
    reference_var_explained = pca$var_explained[seq_along(mediators)],
    point_estimates = point_estimates
  )
}

generate_permutations <- function(k) {
  if (k == 1) {
    return(matrix(1, nrow = 1))
  }

  do.call(
    rbind,
    lapply(seq_len(k), function(first) {
      rest <- generate_permutations(k - 1)
      remaining <- setdiff(seq_len(k), first)
      cbind(first, matrix(remaining[rest], nrow = nrow(rest)))
    })
  )
}

match_pca_components <- function(reference_loadings,
                                 bootstrap_loadings,
                                 bootstrap_scores,
                                 bootstrap_var_explained,
                                 mediator_names = c("PC1_R", "PC2_R", "PC3"),
                                 bootstrap_id = NA_integer_) {
  common_genes <- intersect(rownames(reference_loadings), rownames(bootstrap_loadings))
  if (length(common_genes) < length(mediator_names)) {
    stop("Too few common genes for PCA loading alignment.", call. = FALSE)
  }

  k <- length(mediator_names)
  ref <- reference_loadings[common_genes, mediator_names, drop = FALSE]
  boot <- bootstrap_loadings[common_genes, seq_len(k), drop = FALSE]

  similarity <- stats::cor(ref, boot)
  permutations <- generate_permutations(k)
  scores <- apply(permutations, 1, function(perm) {
    sum(abs(similarity[cbind(seq_len(k), perm)]))
  })
  best_perm <- permutations[which.max(scores), ]

  aligned_loadings <- matrix(NA_real_, nrow = nrow(bootstrap_loadings), ncol = k)
  aligned_scores <- matrix(NA_real_, nrow = nrow(bootstrap_scores), ncol = k)
  colnames(aligned_loadings) <- mediator_names
  colnames(aligned_scores) <- mediator_names
  rownames(aligned_loadings) <- rownames(bootstrap_loadings)
  rownames(aligned_scores) <- rownames(bootstrap_scores)

  diagnostic_rows <- vector("list", k)
  for (j in seq_len(k)) {
    boot_component <- best_perm[[j]]
    loading_correlation <- similarity[j, boot_component]
    sign_multiplier <- ifelse(loading_correlation < 0, -1, 1)

    aligned_loadings[, j] <- sign_multiplier * bootstrap_loadings[, boot_component]
    aligned_scores[, j] <- sign_multiplier * bootstrap_scores[, boot_component]

    diagnostic_rows[[j]] <- tibble::tibble(
      bootstrap_id = bootstrap_id,
      reference_component = mediator_names[[j]],
      matched_bootstrap_component = paste0("PC", boot_component),
      loading_correlation_before_sign = loading_correlation,
      absolute_loading_correlation = abs(loading_correlation),
      sign_flipped = loading_correlation < 0,
      explained_variance_bootstrap = bootstrap_var_explained[[boot_component]]
    )
  }
  diagnostics <- dplyr::bind_rows(diagnostic_rows)

  ref_basis <- qr.Q(qr(as.matrix(ref)))
  boot_basis <- qr.Q(qr(as.matrix(boot[, best_perm, drop = FALSE])))
  singular_values <- svd(t(ref_basis) %*% boot_basis, nu = 0, nv = 0)$d
  singular_values <- pmin(1, pmax(0, singular_values))
  principal_angles_degrees <- acos(singular_values) * 180 / pi

  list(
    aligned_loadings = aligned_loadings,
    aligned_scores = aligned_scores,
    diagnostics = diagnostics,
    similarity_matrix = similarity,
    permutation = best_perm,
    subspace = tibble::tibble(
      bootstrap_id = bootstrap_id,
      subspace_dimension = seq_along(singular_values),
      singular_value = singular_values,
      principal_angle_degrees = principal_angles_degrees
    )
  )
}

validate_pca_alignment_helper <- function() {
  set.seed(11)
  ref <- qr.Q(qr(matrix(stats::rnorm(90), nrow = 30, ncol = 3)))
  colnames(ref) <- c("PC1_R", "PC2_R", "PC3")
  rownames(ref) <- paste0("gene_", seq_len(nrow(ref)))

  perm <- c(3, 1, 2)
  signs <- c(-1, 1, -1)
  boot <- ref[, perm] %*% diag(signs)
  boot <- boot + matrix(stats::rnorm(length(boot), sd = 1e-8), nrow = nrow(boot))
  colnames(boot) <- paste0("PC", 1:3)
  rownames(boot) <- rownames(ref)
  scores <- matrix(stats::rnorm(45), nrow = 15, ncol = 3)

  aligned <- match_pca_components(
    reference_loadings = ref,
    bootstrap_loadings = boot,
    bootstrap_scores = scores,
    bootstrap_var_explained = c(0.4, 0.3, 0.2),
    bootstrap_id = 0
  )

  list(
    passed = all(aligned$diagnostics$absolute_loading_correlation > 0.999),
    diagnostics = aligned$diagnostics,
    similarity_matrix = aligned$similarity_matrix
  )
}

summarize_refit_bootstrap_quantities <- function(decomposition,
                                                bootstrap_id,
                                                high_quality_alignment = TRUE) {
  summary_draws <- decomposition$summary |>
    dplyr::mutate(
      bootstrap_id = bootstrap_id,
      contrast = pca_refit_contrast_label(exposure),
      direct_effect = NDE,
      total_NIE = NIE,
      high_quality_alignment = high_quality_alignment
    ) |>
    dplyr::select(
      bootstrap_id,
      high_quality_alignment,
      contrast,
      exposure,
      TE,
      direct_effect,
      total_NIE,
      PM
    )

  path_draws <- decomposition$path |>
    dplyr::mutate(
      bootstrap_id = bootstrap_id,
      contrast = pca_refit_contrast_label(exposure),
      alpha = alpha_X_to_M,
      beta = beta_M_to_Y,
      high_quality_alignment = high_quality_alignment
    ) |>
    dplyr::select(
      bootstrap_id,
      high_quality_alignment,
      contrast,
      exposure,
      mediator,
      alpha,
      beta,
      indirect_component
    )

  quantity_draws <- dplyr::bind_rows(
    summary_draws |>
      tidyr::pivot_longer(
        cols = c(TE, direct_effect, total_NIE, PM),
        names_to = "quantity",
        values_to = "value"
      ) |>
      dplyr::mutate(mediator = NA_character_),
    path_draws |>
      dplyr::transmute(
        bootstrap_id,
        high_quality_alignment,
        contrast,
        exposure,
        mediator,
        quantity = paste0("IE_", mediator),
        value = indirect_component
      )
  ) |>
    dplyr::select(bootstrap_id, high_quality_alignment, contrast, exposure, mediator, quantity, value)

  list(
    summary_draws = summary_draws,
    path_draws = path_draws,
    quantity_draws = quantity_draws
  )
}

run_one_pca_refit_bootstrap <- function(bootstrap_id,
                                        replicate_seed,
                                        reference,
                                        min_loading_correlation = 0.70,
                                        mediators = c("PC1_R", "PC2_R", "PC3")) {
  start_time <- proc.time()

  tryCatch({
    set.seed(replicate_seed)
    n <- nrow(reference$dat_full_allpc)
    idx <- sample(seq_len(n), size = n, replace = TRUE)

    dat_boot_meta <- reference$dat_full_allpc[idx, c("Y", "X_cat", "X_adjacent", "X_far", "x_coord", "y_coord")]
    M_expr_boot <- reference$M_expr_analysis[idx, , drop = FALSE]

    pca_elapsed <- system.time({
      pca_boot <- fit_pca_mediators(M_expr_boot)
    })

    bootstrap_loadings <- as.matrix(pca_boot$pca_fit$rotation[, seq_along(mediators), drop = FALSE])
    bootstrap_scores <- as.matrix(pca_boot$pca_fit$x[, seq_along(mediators), drop = FALSE])

    alignment <- match_pca_components(
      reference_loadings = reference$reference_loadings,
      bootstrap_loadings = bootstrap_loadings,
      bootstrap_scores = bootstrap_scores,
      bootstrap_var_explained = pca_boot$var_explained,
      mediator_names = mediators,
      bootstrap_id = bootstrap_id
    )

    alignment_warning <- any(alignment$diagnostics$absolute_loading_correlation < min_loading_correlation)
    high_quality_alignment <- !alignment_warning

    dat_boot <- dplyr::bind_cols(
      dat_boot_meta,
      as.data.frame(alignment$aligned_scores[, mediators, drop = FALSE])
    )

    mediation_elapsed <- system.time({
      decomposition <- decompose_linear_multix(
        data = dat_boot,
        exposures = c("X_adjacent", "X_far"),
        outcome = "Y",
        mediators = mediators,
        covariates = c("x_coord", "y_coord")
      )
    })

    if (any(!is.finite(decomposition$summary$TE)) ||
        any(!is.finite(decomposition$summary$NDE)) ||
        any(!is.finite(decomposition$summary$NIE)) ||
        any(!is.finite(decomposition$path$alpha_X_to_M)) ||
        any(!is.finite(decomposition$path$beta_M_to_Y)) ||
        any(!is.finite(decomposition$path$indirect_component))) {
      stop("Non-finite mediation coefficient or decomposition quantity.", call. = FALSE)
    }

    draws <- summarize_refit_bootstrap_quantities(
      decomposition = decomposition,
      bootstrap_id = bootstrap_id,
      high_quality_alignment = high_quality_alignment
    )

    elapsed <- proc.time() - start_time
    list(
      success = TRUE,
      bootstrap_id = bootstrap_id,
      failure_reason = NA_character_,
      summary_draws = draws$summary_draws,
      path_draws = draws$path_draws,
      quantity_draws = draws$quantity_draws,
      alignment_diagnostics = alignment$diagnostics |>
        dplyr::mutate(alignment_warning = alignment_warning),
      alignment_summary = tibble::tibble(
        bootstrap_id = bootstrap_id,
        min_matched_abs_correlation = min(alignment$diagnostics$absolute_loading_correlation),
        mean_matched_abs_correlation = mean(alignment$diagnostics$absolute_loading_correlation),
        alignment_warning = alignment_warning,
        high_quality_alignment = high_quality_alignment
      ),
      subspace_diagnostics = alignment$subspace,
      timing = tibble::tibble(
        bootstrap_id = bootstrap_id,
        pca_elapsed_seconds = unname(pca_elapsed[["elapsed"]]),
        mediation_elapsed_seconds = unname(mediation_elapsed[["elapsed"]]),
        total_elapsed_seconds = unname(elapsed[["elapsed"]])
      )
    )
  }, error = function(e) {
    elapsed <- proc.time() - start_time
    list(
      success = FALSE,
      bootstrap_id = bootstrap_id,
      failure_reason = conditionMessage(e),
      summary_draws = tibble::tibble(),
      path_draws = tibble::tibble(),
      quantity_draws = tibble::tibble(),
      alignment_diagnostics = tibble::tibble(),
      alignment_summary = tibble::tibble(
        bootstrap_id = bootstrap_id,
        min_matched_abs_correlation = NA_real_,
        mean_matched_abs_correlation = NA_real_,
        alignment_warning = NA,
        high_quality_alignment = FALSE
      ),
      subspace_diagnostics = tibble::tibble(),
      timing = tibble::tibble(
        bootstrap_id = bootstrap_id,
        pca_elapsed_seconds = NA_real_,
        mediation_elapsed_seconds = NA_real_,
        total_elapsed_seconds = unname(elapsed[["elapsed"]])
      )
    )
  })
}

run_pca_refit_bootstrap_replicates <- function(reference,
                                               B = 1000,
                                               seed = 123,
                                               n_cores = 1,
                                               min_loading_correlation = 0.70,
                                               mediators = c("PC1_R", "PC2_R", "PC3")) {
  set.seed(seed)
  replicate_seeds <- sample.int(.Machine$integer.max, B)

  run_one <- function(b) {
    run_one_pca_refit_bootstrap(
      bootstrap_id = b,
      replicate_seed = replicate_seeds[[b]],
      reference = reference,
      min_loading_correlation = min_loading_correlation,
      mediators = mediators
    )
  }

  if (n_cores > 1 && .Platform$OS.type != "windows") {
    results <- parallel::mclapply(seq_len(B), run_one, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    if (n_cores > 1 && .Platform$OS.type == "windows") {
      message("Parallel mclapply is unavailable on Windows; using serial mode.")
    }
    results <- lapply(seq_len(B), run_one)
  }

  results
}

empty_bootstrap_table <- function(table_name) {
  switch(
    table_name,
    summary_draws = tibble::tibble(
      bootstrap_id = integer(),
      high_quality_alignment = logical(),
      contrast = character(),
      exposure = character(),
      TE = double(),
      direct_effect = double(),
      total_NIE = double(),
      PM = double()
    ),
    path_draws = tibble::tibble(
      bootstrap_id = integer(),
      high_quality_alignment = logical(),
      contrast = character(),
      exposure = character(),
      mediator = character(),
      alpha = double(),
      beta = double(),
      indirect_component = double()
    ),
    quantity_draws = tibble::tibble(
      bootstrap_id = integer(),
      high_quality_alignment = logical(),
      contrast = character(),
      exposure = character(),
      mediator = character(),
      quantity = character(),
      value = double()
    ),
    alignment_diagnostics = tibble::tibble(
      bootstrap_id = integer(),
      reference_component = character(),
      matched_bootstrap_component = character(),
      loading_correlation_before_sign = double(),
      absolute_loading_correlation = double(),
      sign_flipped = logical(),
      explained_variance_bootstrap = double(),
      alignment_warning = logical()
    ),
    alignment_summary = tibble::tibble(
      bootstrap_id = integer(),
      min_matched_abs_correlation = double(),
      mean_matched_abs_correlation = double(),
      alignment_warning = logical(),
      high_quality_alignment = logical()
    ),
    subspace_diagnostics = tibble::tibble(
      bootstrap_id = integer(),
      subspace_dimension = integer(),
      singular_value = double(),
      principal_angle_degrees = double()
    ),
    timing = tibble::tibble(
      bootstrap_id = integer(),
      pca_elapsed_seconds = double(),
      mediation_elapsed_seconds = double(),
      total_elapsed_seconds = double()
    ),
    tibble::tibble()
  )
}

bind_successful_bootstrap_table <- function(replicate_results, table_name) {
  bound <- dplyr::bind_rows(purrr::map(replicate_results, table_name))
  if (nrow(bound) == 0 && ncol(bound) == 0) {
    return(empty_bootstrap_table(table_name))
  }
  bound
}

summarize_quantity_draws <- function(quantity_draws,
                                    value_col_prefix,
                                    high_quality_only = FALSE) {
  data <- quantity_draws
  if (isTRUE(high_quality_only) && "high_quality_alignment" %in% names(data)) {
    data <- dplyr::filter(data, high_quality_alignment)
  }

  if (nrow(data) == 0) {
    return(tibble::tibble(
      contrast = character(),
      exposure = character(),
      mediator = character(),
      quantity = character(),
      "{value_col_prefix}_boot_mean" := double(),
      "{value_col_prefix}_q025" := double(),
      "{value_col_prefix}_q975" := double(),
      "{value_col_prefix}_SD" := double(),
      "{value_col_prefix}_B" := integer()
    ))
  }

  data |>
    dplyr::group_by(contrast, exposure, mediator, quantity) |>
    dplyr::summarise(
      "{value_col_prefix}_boot_mean" := mean(value, na.rm = TRUE),
      "{value_col_prefix}_q025" := stats::quantile(value, 0.025, na.rm = TRUE),
      "{value_col_prefix}_q975" := stats::quantile(value, 0.975, na.rm = TRUE),
      "{value_col_prefix}_SD" := stats::sd(value, na.rm = TRUE),
      "{value_col_prefix}_B" := dplyr::n(),
      .groups = "drop"
    )
}

make_point_quantity_table <- function(point_estimates) {
  summary_quantities <- point_estimates$summary |>
    dplyr::mutate(
      contrast = pca_refit_contrast_label(exposure),
      direct_effect = NDE,
      total_NIE = NIE
    ) |>
    dplyr::select(contrast, exposure, TE, direct_effect, total_NIE, PM) |>
    tidyr::pivot_longer(
      cols = c(TE, direct_effect, total_NIE, PM),
      names_to = "quantity",
      values_to = "point_estimate"
    ) |>
    dplyr::mutate(mediator = NA_character_)

  path_quantities <- point_estimates$path |>
    dplyr::mutate(
      contrast = pca_refit_contrast_label(exposure),
      quantity = paste0("IE_", mediator),
      point_estimate = indirect_component
    ) |>
    dplyr::select(contrast, exposure, mediator, quantity, point_estimate)

  dplyr::bind_rows(summary_quantities, path_quantities) |>
    dplyr::select(contrast, exposure, mediator, quantity, point_estimate)
}

make_fixed_pca_quantity_draws <- function(fixed_bootstrap) {
  summary_draws <- fixed_bootstrap$boot_summary_df |>
    dplyr::mutate(
      bootstrap_id = as.integer(boot_id),
      contrast = pca_refit_contrast_label(exposure),
      direct_effect = NDE,
      total_NIE = NIE
    ) |>
    dplyr::select(bootstrap_id, contrast, exposure, TE, direct_effect, total_NIE, PM) |>
    tidyr::pivot_longer(
      cols = c(TE, direct_effect, total_NIE, PM),
      names_to = "quantity",
      values_to = "value"
    ) |>
    dplyr::mutate(mediator = NA_character_)

  path_draws <- fixed_bootstrap$boot_path_df |>
    dplyr::mutate(
      bootstrap_id = as.integer(boot_id),
      contrast = pca_refit_contrast_label(exposure),
      quantity = paste0("IE_", mediator),
      value = indirect_component
    ) |>
    dplyr::select(bootstrap_id, contrast, exposure, mediator, quantity, value)

  dplyr::bind_rows(summary_draws, path_draws) |>
    dplyr::select(bootstrap_id, contrast, exposure, mediator, quantity, value)
}

compare_fixed_vs_refit_bootstrap <- function(point_estimates,
                                             fixed_quantity_draws,
                                             refit_quantity_draws,
                                             high_quality_only = FALSE) {
  point_tbl <- make_point_quantity_table(point_estimates)
  fixed_summary <- summarize_quantity_draws(
    quantity_draws = fixed_quantity_draws,
    value_col_prefix = "fixed_PCA"
  )
  refit_summary <- summarize_quantity_draws(
    quantity_draws = refit_quantity_draws,
    value_col_prefix = "refit_PCA",
    high_quality_only = high_quality_only
  )

  point_tbl |>
    dplyr::left_join(fixed_summary, by = c("contrast", "exposure", "mediator", "quantity")) |>
    dplyr::left_join(refit_summary, by = c("contrast", "exposure", "mediator", "quantity")) |>
    dplyr::mutate(
      CI_width_fixed = fixed_PCA_q975 - fixed_PCA_q025,
      CI_width_refit = refit_PCA_q975 - refit_PCA_q025,
      CI_width_ratio = CI_width_refit / CI_width_fixed,
      refit_summary = ifelse(high_quality_only, "high_quality_only", "all_successful_replicates")
    ) |>
    dplyr::select(
      contrast,
      mediator,
      quantity,
      point_estimate,
      fixed_PCA_boot_mean,
      fixed_PCA_q025,
      fixed_PCA_q975,
      fixed_PCA_SD,
      refit_PCA_boot_mean,
      refit_PCA_q025,
      refit_PCA_q975,
      refit_PCA_SD,
      CI_width_fixed,
      CI_width_refit,
      CI_width_ratio,
      fixed_PCA_B,
      refit_PCA_B,
      refit_summary
    )
}

summarize_explained_variance_stability <- function(alignment_diagnostics) {
  if (nrow(alignment_diagnostics) == 0) {
    return(tibble::tibble(
      reference_component = character(),
      mean = double(),
      SD = double(),
      q025 = double(),
      median = double(),
      q975 = double()
    ))
  }

  alignment_diagnostics |>
    dplyr::group_by(reference_component) |>
    dplyr::summarise(
      mean = mean(explained_variance_bootstrap, na.rm = TRUE),
      SD = stats::sd(explained_variance_bootstrap, na.rm = TRUE),
      q025 = stats::quantile(explained_variance_bootstrap, 0.025, na.rm = TRUE),
      median = stats::median(explained_variance_bootstrap, na.rm = TRUE),
      q975 = stats::quantile(explained_variance_bootstrap, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_loading_alignment_stability <- function(alignment_diagnostics) {
  if (nrow(alignment_diagnostics) == 0) {
    return(tibble::tibble(
      reference_component = character(),
      mean_abs_loading_correlation = double(),
      SD_abs_loading_correlation = double(),
      min_abs_loading_correlation = double(),
      q025_abs_loading_correlation = double(),
      median_abs_loading_correlation = double(),
      q975_abs_loading_correlation = double(),
      sign_flip_rate = double()
    ))
  }

  alignment_diagnostics |>
    dplyr::group_by(reference_component) |>
    dplyr::summarise(
      mean_abs_loading_correlation = mean(absolute_loading_correlation, na.rm = TRUE),
      SD_abs_loading_correlation = stats::sd(absolute_loading_correlation, na.rm = TRUE),
      min_abs_loading_correlation = min(absolute_loading_correlation, na.rm = TRUE),
      q025_abs_loading_correlation = stats::quantile(absolute_loading_correlation, 0.025, na.rm = TRUE),
      median_abs_loading_correlation = stats::median(absolute_loading_correlation, na.rm = TRUE),
      q975_abs_loading_correlation = stats::quantile(absolute_loading_correlation, 0.975, na.rm = TRUE),
      sign_flip_rate = mean(sign_flipped, na.rm = TRUE),
      .groups = "drop"
    )
}

make_failure_summary <- function(replicate_results,
                                 requested_B,
                                 alignment_summary) {
  success <- vapply(replicate_results, function(x) isTRUE(x$success), logical(1))
  failures <- tibble::tibble(
    bootstrap_id = vapply(replicate_results, function(x) x$bootstrap_id, integer(1)),
    success = success,
    failure_reason = vapply(replicate_results, function(x) x$failure_reason %||% NA_character_, character(1))
  )

  tibble::tibble(
    requested_B = requested_B,
    successful_B = sum(success),
    failed_B = sum(!success),
    alignment_warning_B = sum(alignment_summary$alignment_warning, na.rm = TRUE)
  ) |>
    dplyr::mutate(failure_details = list(failures))
}

run_pca_refit_bootstrap_analysis <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    B = 1000,
    seed = 123,
    n_cores = 1,
    min_loading_correlation = 0.70,
    run_fixed_pca_bootstrap = TRUE) {
  check_pca_refit_bootstrap_packages()

  total_elapsed <- system.time({
    reference <- prepare_pca_refit_bootstrap_reference(rds_path = rds_path)

    alignment_sanity_check <- validate_pca_alignment_helper()

    fixed_bootstrap <- NULL
    fixed_quantity_draws <- tibble::tibble()
    fixed_elapsed <- NA_real_
    if (isTRUE(run_fixed_pca_bootstrap)) {
      fixed_time <- system.time({
        fixed_bootstrap <- bootstrap_decomposition(
          data = reference$dat_full_allpc,
          B = B,
          seed = seed,
          exposures = c("X_adjacent", "X_far"),
          outcome = "Y",
          mediators = c("PC1_R", "PC2_R", "PC3"),
          covariates = c("x_coord", "y_coord")
        )
        fixed_quantity_draws <- make_fixed_pca_quantity_draws(fixed_bootstrap)
      })
      fixed_elapsed <- unname(fixed_time[["elapsed"]])
    }

    refit_time <- system.time({
      replicate_results <- run_pca_refit_bootstrap_replicates(
        reference = reference,
        B = B,
        seed = seed,
        n_cores = n_cores,
        min_loading_correlation = min_loading_correlation
      )
    })

    summary_draws <- bind_successful_bootstrap_table(replicate_results, "summary_draws")
    path_draws <- bind_successful_bootstrap_table(replicate_results, "path_draws")
    quantity_draws <- bind_successful_bootstrap_table(replicate_results, "quantity_draws")
    alignment_diagnostics <- bind_successful_bootstrap_table(replicate_results, "alignment_diagnostics")
    alignment_summary <- bind_successful_bootstrap_table(replicate_results, "alignment_summary")
    subspace_diagnostics <- bind_successful_bootstrap_table(replicate_results, "subspace_diagnostics")
    replicate_timing <- bind_successful_bootstrap_table(replicate_results, "timing")

    failure_summary <- make_failure_summary(
      replicate_results = replicate_results,
      requested_B = B,
      alignment_summary = alignment_summary
    )

    refit_quantity_summary_all <- summarize_quantity_draws(
      quantity_draws = quantity_draws,
      value_col_prefix = "refit_PCA",
      high_quality_only = FALSE
    )
    refit_quantity_summary_high_quality <- summarize_quantity_draws(
      quantity_draws = quantity_draws,
      value_col_prefix = "refit_PCA",
      high_quality_only = TRUE
    )

    fixed_vs_refit_comparison <- NULL
    fixed_vs_refit_comparison_high_quality <- NULL
    if (isTRUE(run_fixed_pca_bootstrap)) {
      fixed_vs_refit_comparison <- compare_fixed_vs_refit_bootstrap(
        point_estimates = reference$point_estimates,
        fixed_quantity_draws = fixed_quantity_draws,
        refit_quantity_draws = quantity_draws,
        high_quality_only = FALSE
      )
      fixed_vs_refit_comparison_high_quality <- compare_fixed_vs_refit_bootstrap(
        point_estimates = reference$point_estimates,
        fixed_quantity_draws = fixed_quantity_draws,
        refit_quantity_draws = quantity_draws,
        high_quality_only = TRUE
      )
    }

    explained_variance_stability <- summarize_explained_variance_stability(alignment_diagnostics)
    loading_alignment_stability <- summarize_loading_alignment_stability(alignment_diagnostics)
  })

  timing <- tibble::tibble(
    component = c(
      "fixed_PCA_bootstrap",
      "refit_PCA_bootstrap",
      "total",
      "average_per_refit_bootstrap",
      "average_pca_fit_per_refit_bootstrap",
      "average_mediation_fit_per_refit_bootstrap"
    ),
    elapsed_seconds = c(
      fixed_elapsed,
      unname(refit_time[["elapsed"]]),
      unname(total_elapsed[["elapsed"]]),
      mean(replicate_timing$total_elapsed_seconds, na.rm = TRUE),
      mean(replicate_timing$pca_elapsed_seconds, na.rm = TRUE),
      mean(replicate_timing$mediation_elapsed_seconds, na.rm = TRUE)
    )
  )

  list(
    config = list(
      rds_path = rds_path,
      B = B,
      seed = seed,
      n_cores = n_cores,
      min_loading_correlation = min_loading_correlation,
      run_fixed_pca_bootstrap = run_fixed_pca_bootstrap,
      resampling_unit = "spot",
      resampling = "sample n analysis spots with replacement"
    ),
    reference = list(
      pca_fit = reference$pca$pca_fit,
      reference_loadings = reference$reference_loadings,
      reference_var_explained = reference$reference_var_explained,
      point_estimates = reference$point_estimates,
      dat_full_allpc = reference$dat_full_allpc
    ),
    validation = list(
      alignment_sanity_check = alignment_sanity_check
    ),
    fixed_pca_bootstrap = fixed_bootstrap,
    bootstrap_draws = list(
      summary = summary_draws,
      path = path_draws,
      quantities = quantity_draws,
      fixed_pca_quantities = fixed_quantity_draws
    ),
    alignment_diagnostics = alignment_diagnostics,
    alignment_summary = alignment_summary,
    subspace_diagnostics = subspace_diagnostics,
    explained_variance_stability = explained_variance_stability,
    loading_alignment_stability = loading_alignment_stability,
    failure_summary = failure_summary,
    summaries = list(
      refit_quantity_summary_all = refit_quantity_summary_all,
      refit_quantity_summary_high_quality = refit_quantity_summary_high_quality
    ),
    fixed_vs_refit_comparison = fixed_vs_refit_comparison,
    fixed_vs_refit_comparison_high_quality = fixed_vs_refit_comparison_high_quality,
    timing = timing
  )
}
