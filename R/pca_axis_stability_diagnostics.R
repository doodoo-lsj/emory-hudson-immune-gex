# Axis stability diagnostics for saved PCA-refit bootstrap results.
#
# The saved PCA-refit object already contains matched loading correlations,
# sign/permutation alignment decisions, explained variance for matched axes, and
# K=3 subspace principal-angle diagnostics. Full 3x3 raw loading similarity
# matrices are only available if the PCA refit is recomputed from the saved
# bootstrap seed/config, because the original bootstrap result intentionally did
# not store all loading matrices for memory reasons.

check_pca_axis_stability_packages <- function() {
  required <- c("dplyr", "tidyr", "tibble", "purrr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "PCA axis stability diagnostics require installed R packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

load_saved_pca_refit_results <- function(path = "pca_refit_bootstrap.RData",
                                         object_name = "pca_refit_bootstrap_results") {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop("Object not found in ", path, ": ", object_name, call. = FALSE)
  }
  get(object_name, envir = env, inherits = FALSE)
}

expected_bootstrap_component <- function(reference_component) {
  dplyr::case_when(
    reference_component == "PC1_R" ~ "PC1",
    reference_component == "PC2_R" ~ "PC2",
    reference_component == "PC3" ~ "PC3",
    TRUE ~ reference_component
  )
}

prepare_axis_alignment_table <- function(alignment_diagnostics) {
  alignment_diagnostics |>
    dplyr::mutate(
      expected_bootstrap_component = expected_bootstrap_component(reference_component),
      permutation_changed = matched_bootstrap_component != expected_bootstrap_component,
      aligned_loading_correlation = abs(loading_correlation_before_sign),
      aligned_abs_loading_correlation = absolute_loading_correlation,
      aligned_abs_loading_cosine = NA_real_,
      note = "Cosine similarity is NA unless full similarity matrices are recomputed."
    )
}

summarize_axis_similarity_distribution <- function(axis_alignment) {
  axis_alignment |>
    dplyr::group_by(reference_component) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean_abs_loading_correlation = mean(aligned_abs_loading_correlation, na.rm = TRUE),
      sd_abs_loading_correlation = stats::sd(aligned_abs_loading_correlation, na.rm = TRUE),
      min_abs_loading_correlation = min(aligned_abs_loading_correlation, na.rm = TRUE),
      q025_abs_loading_correlation = stats::quantile(aligned_abs_loading_correlation, 0.025, na.rm = TRUE),
      median_abs_loading_correlation = stats::median(aligned_abs_loading_correlation, na.rm = TRUE),
      q975_abs_loading_correlation = stats::quantile(aligned_abs_loading_correlation, 0.975, na.rm = TRUE),
      max_abs_loading_correlation = max(aligned_abs_loading_correlation, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_permutation_events <- function(axis_alignment) {
  per_axis <- axis_alignment |>
    dplyr::count(
      reference_component,
      expected_bootstrap_component,
      matched_bootstrap_component,
      permutation_changed,
      name = "n"
    ) |>
    dplyr::group_by(reference_component) |>
    dplyr::mutate(proportion = n / sum(n)) |>
    dplyr::ungroup()

  per_replicate <- axis_alignment |>
    dplyr::arrange(bootstrap_id, reference_component) |>
    dplyr::group_by(bootstrap_id) |>
    dplyr::summarise(
      permutation_pattern = paste(matched_bootstrap_component, collapse = " -> "),
      any_permutation_changed = any(permutation_changed),
      n_axes_permuted = sum(permutation_changed),
      .groups = "drop"
    )

  list(
    per_axis = per_axis,
    per_replicate = per_replicate,
    total = per_replicate |>
      dplyr::summarise(
        bootstrap_replicates = dplyr::n(),
        swap_replicates = sum(any_permutation_changed),
        no_swap_replicates = sum(!any_permutation_changed),
        swap_rate = mean(any_permutation_changed)
      ),
    pattern_counts = per_replicate |>
      dplyr::count(permutation_pattern, any_permutation_changed, name = "n") |>
      dplyr::mutate(proportion = n / sum(n))
  )
}

summarize_sign_flips <- function(axis_alignment) {
  per_axis <- axis_alignment |>
    dplyr::group_by(reference_component) |>
    dplyr::summarise(
      n = dplyr::n(),
      sign_flip_n = sum(sign_flipped, na.rm = TRUE),
      sign_flip_rate = mean(sign_flipped, na.rm = TRUE),
      .groups = "drop"
    )

  per_replicate <- axis_alignment |>
    dplyr::group_by(bootstrap_id) |>
    dplyr::summarise(
      any_sign_flip = any(sign_flipped, na.rm = TRUE),
      n_sign_flipped_axes = sum(sign_flipped, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    per_axis = per_axis,
    per_replicate = per_replicate,
    total = per_replicate |>
      dplyr::summarise(
        bootstrap_replicates = dplyr::n(),
        sign_flip_replicates = sum(any_sign_flip),
        sign_flip_replicate_rate = mean(any_sign_flip),
        total_axis_sign_flips = sum(n_sign_flipped_axes)
      )
  )
}

summarize_explained_variance_ordering <- function(axis_alignment) {
  ev_wide <- axis_alignment |>
    dplyr::select(bootstrap_id, reference_component, explained_variance_bootstrap) |>
    tidyr::pivot_wider(
      names_from = reference_component,
      values_from = explained_variance_bootstrap
    ) |>
    dplyr::mutate(
      gap_PC1_R_minus_PC2_R = PC1_R - PC2_R,
      gap_PC2_R_minus_PC3 = PC2_R - PC3,
      gap_PC1_R_minus_PC3 = PC1_R - PC3,
      aligned_ev_order_preserved = PC1_R >= PC2_R & PC2_R >= PC3
    )

  component_summary <- axis_alignment |>
    dplyr::group_by(reference_component) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean_explained_variance = mean(explained_variance_bootstrap, na.rm = TRUE),
      sd_explained_variance = stats::sd(explained_variance_bootstrap, na.rm = TRUE),
      q025_explained_variance = stats::quantile(explained_variance_bootstrap, 0.025, na.rm = TRUE),
      median_explained_variance = stats::median(explained_variance_bootstrap, na.rm = TRUE),
      q975_explained_variance = stats::quantile(explained_variance_bootstrap, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  gap_summary <- ev_wide |>
    tidyr::pivot_longer(
      cols = dplyr::starts_with("gap_"),
      names_to = "gap",
      values_to = "value"
    ) |>
    dplyr::group_by(gap) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(value, na.rm = TRUE),
      sd = stats::sd(value, na.rm = TRUE),
      q025 = stats::quantile(value, 0.025, na.rm = TRUE),
      median = stats::median(value, na.rm = TRUE),
      q975 = stats::quantile(value, 0.975, na.rm = TRUE),
      proportion_positive = mean(value > 0, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    per_replicate = ev_wide,
    component_summary = component_summary,
    gap_summary = gap_summary,
    ordering_summary = ev_wide |>
      dplyr::summarise(
        bootstrap_replicates = dplyr::n(),
        aligned_ev_order_preserved_n = sum(aligned_ev_order_preserved),
        aligned_ev_order_preserved_rate = mean(aligned_ev_order_preserved)
      )
  )
}

summarize_subspace_stability <- function(subspace_diagnostics) {
  per_dimension <- subspace_diagnostics |>
    dplyr::group_by(subspace_dimension) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean_singular_value = mean(singular_value, na.rm = TRUE),
      sd_singular_value = stats::sd(singular_value, na.rm = TRUE),
      q025_singular_value = stats::quantile(singular_value, 0.025, na.rm = TRUE),
      median_singular_value = stats::median(singular_value, na.rm = TRUE),
      q975_singular_value = stats::quantile(singular_value, 0.975, na.rm = TRUE),
      mean_principal_angle_degrees = mean(principal_angle_degrees, na.rm = TRUE),
      q975_principal_angle_degrees = stats::quantile(principal_angle_degrees, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  per_replicate <- subspace_diagnostics |>
    dplyr::group_by(bootstrap_id) |>
    dplyr::summarise(
      min_subspace_singular_value = min(singular_value, na.rm = TRUE),
      mean_subspace_singular_value = mean(singular_value, na.rm = TRUE),
      max_principal_angle_degrees = max(principal_angle_degrees, na.rm = TRUE),
      mean_principal_angle_degrees = mean(principal_angle_degrees, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    per_dimension = per_dimension,
    per_replicate = per_replicate,
    replicate_summary = per_replicate |>
      dplyr::summarise(
        n = dplyr::n(),
        mean_min_subspace_singular_value = mean(min_subspace_singular_value, na.rm = TRUE),
        q025_min_subspace_singular_value = stats::quantile(min_subspace_singular_value, 0.025, na.rm = TRUE),
        median_min_subspace_singular_value = stats::median(min_subspace_singular_value, na.rm = TRUE),
        mean_max_principal_angle_degrees = mean(max_principal_angle_degrees, na.rm = TRUE),
        q975_max_principal_angle_degrees = stats::quantile(max_principal_angle_degrees, 0.975, na.rm = TRUE)
      )
  )
}

compute_axis_vs_subspace_summary <- function(axis_similarity_summary,
                                             subspace_summary) {
  worst_axis <- axis_similarity_summary |>
    dplyr::arrange(median_abs_loading_correlation) |>
    dplyr::slice(1)

  tibble::tibble(
    lowest_median_axis_component = worst_axis$reference_component,
    lowest_median_axis_abs_loading_correlation = worst_axis$median_abs_loading_correlation,
    median_min_subspace_singular_value = subspace_summary$replicate_summary$median_min_subspace_singular_value,
    q025_min_subspace_singular_value = subspace_summary$replicate_summary$q025_min_subspace_singular_value,
    interpretation_guardrail = paste(
      "Compare individual-axis similarity with K=3 subspace singular values.",
      "Low individual PC similarity with high subspace similarity indicates axis rotation within a stable 3D representation."
    )
  )
}

unit_normalize_columns <- function(x) {
  norms <- sqrt(colSums(x^2))
  sweep(x, 2, norms, "/")
}

recompute_similarity_matrix_for_bootstrap <- function(bootstrap_id,
                                                      replicate_seed,
                                                      reference,
                                                      mediators = c("PC1_R", "PC2_R", "PC3")) {
  set.seed(replicate_seed)
  n <- nrow(reference$dat_full_allpc)
  idx <- sample(seq_len(n), size = n, replace = TRUE)
  M_expr_boot <- reference$M_expr_analysis[idx, , drop = FALSE]
  pca_boot <- fit_pca_mediators(M_expr_boot)

  bootstrap_loadings <- as.matrix(pca_boot$pca_fit$rotation[, seq_along(mediators), drop = FALSE])
  colnames(bootstrap_loadings) <- paste0("PC", seq_along(mediators))

  common_genes <- intersect(rownames(reference$reference_loadings), rownames(bootstrap_loadings))
  ref <- reference$reference_loadings[common_genes, mediators, drop = FALSE]
  boot <- bootstrap_loadings[common_genes, , drop = FALSE]

  loading_correlation <- stats::cor(ref, boot)
  loading_cosine <- t(unit_normalize_columns(ref)) %*% unit_normalize_columns(boot)

  correlation_long <- as.data.frame(loading_correlation) |>
    tibble::rownames_to_column("reference_component") |>
    tidyr::pivot_longer(
      cols = -reference_component,
      names_to = "bootstrap_component",
      values_to = "loading_correlation"
    )

  cosine_long <- as.data.frame(loading_cosine) |>
    tibble::rownames_to_column("reference_component") |>
    tidyr::pivot_longer(
      cols = -reference_component,
      names_to = "bootstrap_component",
      values_to = "loading_cosine"
    )

  dplyr::left_join(correlation_long, cosine_long, by = c("reference_component", "bootstrap_component")) |>
    dplyr::mutate(
      bootstrap_id = bootstrap_id,
      absolute_loading_correlation = abs(loading_correlation),
      absolute_loading_cosine = abs(loading_cosine),
      .before = 1
    )
}

recompute_similarity_matrices <- function(results,
                                          bootstrap_ids = NULL,
                                          n_cores = 1,
                                          mediators = c("PC1_R", "PC2_R", "PC3")) {
  if (!exists("prepare_pca_refit_bootstrap_reference", mode = "function") ||
      !exists("fit_pca_mediators", mode = "function")) {
    stop("Source R/pca_mediation_pipeline.R and R/pca_refit_bootstrap.R before recomputing similarity matrices.", call. = FALSE)
  }

  reference <- prepare_pca_refit_bootstrap_reference(
    rds_path = results$config$rds_path,
    mediators = mediators
  )
  set.seed(results$config$seed)
  replicate_seeds <- sample.int(.Machine$integer.max, results$config$B)

  if (is.null(bootstrap_ids)) {
    bootstrap_ids <- sort(unique(results$alignment_diagnostics$bootstrap_id))
  }

  run_one <- function(b) {
    recompute_similarity_matrix_for_bootstrap(
      bootstrap_id = b,
      replicate_seed = replicate_seeds[[b]],
      reference = reference,
      mediators = mediators
    )
  }

  if (n_cores > 1 && .Platform$OS.type != "windows") {
    dplyr::bind_rows(parallel::mclapply(bootstrap_ids, run_one, mc.cores = n_cores, mc.preschedule = FALSE))
  } else {
    dplyr::bind_rows(lapply(bootstrap_ids, run_one))
  }
}

run_pca_axis_stability_diagnostics <- function(results = NULL,
                                               results_path = "pca_refit_bootstrap.RData",
                                               object_name = "pca_refit_bootstrap_results",
                                               recompute_full_similarity_matrices = FALSE,
                                               bootstrap_ids_for_recompute = NULL,
                                               n_cores = 1) {
  check_pca_axis_stability_packages()

  if (is.null(results)) {
    results <- load_saved_pca_refit_results(path = results_path, object_name = object_name)
  }

  axis_alignment <- prepare_axis_alignment_table(results$alignment_diagnostics)
  axis_similarity_distribution <- summarize_axis_similarity_distribution(axis_alignment)
  permutation_summary <- summarize_permutation_events(axis_alignment)
  sign_flip_summary <- summarize_sign_flips(axis_alignment)
  explained_variance_ordering <- summarize_explained_variance_ordering(axis_alignment)
  subspace_stability <- summarize_subspace_stability(results$subspace_diagnostics)
  axis_vs_subspace_summary <- compute_axis_vs_subspace_summary(
    axis_similarity_summary = axis_similarity_distribution,
    subspace_summary = subspace_stability
  )

  full_similarity_matrices <- NULL
  full_similarity_note <- paste(
    "Full 3x3 loading correlation/cosine matrices were not stored in the saved PCA-refit object.",
    "Set recompute_full_similarity_matrices = TRUE to reconstruct them by refitting PCA from the saved seed/config."
  )

  if (isTRUE(recompute_full_similarity_matrices)) {
    full_similarity_matrices <- recompute_similarity_matrices(
      results = results,
      bootstrap_ids = bootstrap_ids_for_recompute,
      n_cores = n_cores
    )
    full_similarity_note <- "Full 3x3 loading correlation/cosine matrices were recomputed from saved seed/config."
  }

  list(
    config = list(
      source_results_path = results_path,
      source_object_name = object_name,
      requested_B = results$config$B,
      successful_B = length(unique(results$alignment_diagnostics$bootstrap_id)),
      recompute_full_similarity_matrices = recompute_full_similarity_matrices,
      bootstrap_ids_for_recompute = bootstrap_ids_for_recompute,
      n_cores = n_cores
    ),
    axis_alignment = axis_alignment,
    axis_similarity_distribution = axis_similarity_distribution,
    permutation_summary = permutation_summary,
    sign_flip_summary = sign_flip_summary,
    explained_variance_ordering = explained_variance_ordering,
    subspace_stability = subspace_stability,
    axis_vs_subspace_summary = axis_vs_subspace_summary,
    full_similarity_matrices = full_similarity_matrices,
    full_similarity_note = full_similarity_note,
    alignment_method = list(
      matching = "For K=3, exhaustive permutation maximizes sum_j abs(cor(W_ref[, j], W_boot[, pi(j)])).",
      sign_alignment = "After matching, bootstrap loading and score signs are multiplied by -1 when the matched loading correlation is negative.",
      reference_orientation = "Reference loadings are the current final PC1_R, PC2_R, PC3 orientation from build_loading_data()."
    )
  )
}
