# Spatial adjustment specification diagnostics for the current PCA mediation analysis.
#
# This helper does not replace the main frequentist pipeline. It fixes the same
# analysis data and PCA mediators, then compares linear, quadratic, and 2D spline
# spatial adjustment specifications.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

check_spatial_sensitivity_packages <- function(require_spline = TRUE) {
  required <- c("dplyr", "tidyr", "tibble", "purrr")
  if (isTRUE(require_spline)) {
    required <- c(required, "mgcv")
  }

  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Spatial sensitivity diagnostics require installed R packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

spatial_contrast_label <- function(exposure) {
  dplyr::case_when(
    exposure == "X_adjacent" ~ "adjacent vs inside",
    exposure == "X_far" ~ "far vs inside",
    TRUE ~ exposure
  )
}

prepare_spatial_sensitivity_data <- function(rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
                                             coord_x = "imagecol",
                                             coord_y = "imagerow") {
  needed_functions <- c(
    "load_spatial_mediation_rds",
    "prepare_analysis_metadata",
    "fit_pca_mediators",
    "build_pca_score_data",
    "build_full_analysis_data"
  )
  missing_functions <- needed_functions[!vapply(needed_functions, exists, logical(1), mode = "function")]
  if (length(missing_functions) > 0) {
    stop(
      "Missing current pipeline functions: ",
      paste(missing_functions, collapse = ", "),
      ". Source R/pca_mediation_pipeline.R before running spatial sensitivity diagnostics.",
      call. = FALSE
    )
  }

  gex <- load_spatial_mediation_rds(rds_path)
  analysis_df <- prepare_analysis_metadata(gex, coord_x = coord_x, coord_y = coord_y)
  pca <- fit_pca_mediators(gex$M_expr)
  pca_df <- build_pca_score_data(analysis_df, pca$pca_fit)
  dat <- build_full_analysis_data(pca_df)

  x_center <- mean(dat$x_coord)
  y_center <- mean(dat$y_coord)
  x_scale <- stats::sd(dat$x_coord)
  y_scale <- stats::sd(dat$y_coord)
  if (!is.finite(x_scale) || x_scale <= 0 || !is.finite(y_scale) || y_scale <= 0) {
    stop("Coordinate SD is zero or invalid.", call. = FALSE)
  }

  dat <- dplyr::mutate(
    dat,
    x_std = (.data$x_coord - x_center) / x_scale,
    y_std = (.data$y_coord - y_center) / y_scale
  )

  predictor_scale_summary <- purrr::map_dfr(
    c("Y", "PC1_R", "PC2_R", "PC3", "x_coord", "y_coord", "x_std", "y_std"),
    function(v) {
      tibble::tibble(
        variable = v,
        mean = mean(dat[[v]], na.rm = TRUE),
        sd = stats::sd(dat[[v]], na.rm = TRUE),
        min = min(dat[[v]], na.rm = TRUE),
        max = max(dat[[v]], na.rm = TRUE)
      )
    }
  )

  list(
    gex = gex,
    analysis_df = analysis_df,
    pca = pca,
    pca_df = pca_df,
    data = dat,
    coordinate_scale_info = tibble::tibble(
      variable = c("x_coord", "y_coord"),
      standardized_variable = c("x_std", "y_std"),
      center = c(x_center, y_center),
      scale = c(x_scale, y_scale)
    ),
    predictor_scale_summary = predictor_scale_summary
  )
}

make_spatial_specifications <- function(spline_k = 50,
                                        spline_basis = "tp",
                                        include = c("linear", "quadratic", "spline")) {
  include <- match.arg(include, choices = c("linear", "quadratic", "spline"), several.ok = TRUE)

  specs <- list(
    linear = list(
      name = "linear",
      label = "Linear spatial adjustment",
      engine = "lm",
      spatial_terms = c("x_coord", "y_coord"),
      formula_suffix = "x_coord + y_coord",
      description = "x_coord + y_coord"
    ),
    quadratic = list(
      name = "quadratic",
      label = "Quadratic spatial adjustment",
      engine = "lm",
      spatial_terms = c("x_std", "y_std", "I(x_std^2)", "I(y_std^2)", "x_std:y_std"),
      formula_suffix = "x_std + y_std + I(x_std^2) + I(y_std^2) + x_std:y_std",
      description = "x_std + y_std + I(x_std^2) + I(y_std^2) + x_std:y_std"
    ),
    spline = list(
      name = "spline",
      label = "2D spatial spline adjustment",
      engine = "gam",
      spatial_terms = paste0("s(x_std, y_std, bs = '", spline_basis, "', k = ", spline_k, ")"),
      formula_suffix = paste0("s(x_std, y_std, bs = '", spline_basis, "', k = ", spline_k, ")"),
      description = paste0("mgcv::s(x_std, y_std, bs = '", spline_basis, "', k = ", spline_k, ")"),
      spline_k = spline_k,
      spline_basis = spline_basis,
      method = "REML"
    )
  )

  specs[include]
}

fit_spatial_regression <- function(response,
                                   parametric_terms,
                                   spatial_spec,
                                   data) {
  rhs <- paste(c(parametric_terms, spatial_spec$formula_suffix), collapse = " + ")
  formula <- stats::as.formula(paste(response, "~", rhs))

  if (identical(spatial_spec$engine, "gam")) {
    check_spatial_sensitivity_packages(require_spline = TRUE)
    return(mgcv::gam(formula, data = data, method = spatial_spec$method %||% "REML"))
  }

  stats::lm(formula, data = data)
}

extract_parametric_coefficients <- function(fit, terms) {
  if (inherits(fit, "gam")) {
    coef_table <- summary(fit)$p.table
  } else {
    coef_table <- summary(fit)$coefficients
  }

  missing_terms <- setdiff(terms, rownames(coef_table))
  if (length(missing_terms) > 0) {
    stop("Missing requested coefficient terms: ", paste(missing_terms, collapse = ", "), call. = FALSE)
  }

  tibble::tibble(
    term = terms,
    estimate = coef_table[terms, "Estimate"],
    std_error = coef_table[terms, "Std. Error"]
  )
}

compute_residual_matrices <- function(mediator_models,
                                      mediators = c("PC1_R", "PC2_R", "PC3")) {
  residual_matrix <- do.call(
    cbind,
    purrr::map(mediators, function(med) stats::resid(mediator_models[[med]], type = "response"))
  )
  colnames(residual_matrix) <- mediators

  n <- nrow(residual_matrix)
  covariance_mle <- crossprod(residual_matrix) / n
  correlation <- stats::cov2cor(covariance_mle)

  list(
    residuals = residual_matrix,
    covariance = covariance_mle,
    correlation = correlation
  )
}

make_spatial_pairwise_residual_summary <- function(spatial_specification,
                                                   covariance_matrix,
                                                   correlation_matrix) {
  pairs <- utils::combn(colnames(correlation_matrix), 2, simplify = FALSE)

  purrr::map_dfr(pairs, function(pair) {
    tibble::tibble(
      spatial_specification = spatial_specification,
      mediator_1 = pair[[1]],
      mediator_2 = pair[[2]],
      residual_covariance = covariance_matrix[pair[[1]], pair[[2]]],
      residual_correlation = correlation_matrix[pair[[1]], pair[[2]]]
    )
  })
}

fit_one_spatial_specification <- function(data,
                                          spatial_spec,
                                          exposures = c("X_adjacent", "X_far"),
                                          mediators = c("PC1_R", "PC2_R", "PC3"),
                                          outcome = "Y") {
  mediator_models <- stats::setNames(
    purrr::map(mediators, function(med) {
      fit_spatial_regression(
        response = med,
        parametric_terms = exposures,
        spatial_spec = spatial_spec,
        data = data
      )
    }),
    mediators
  )

  outcome_model <- fit_spatial_regression(
    response = outcome,
    parametric_terms = c(exposures, mediators),
    spatial_spec = spatial_spec,
    data = data
  )

  total_model <- fit_spatial_regression(
    response = outcome,
    parametric_terms = exposures,
    spatial_spec = spatial_spec,
    data = data
  )

  residuals <- compute_residual_matrices(mediator_models, mediators = mediators)

  alpha_tbl <- purrr::map_dfr(mediators, function(med) {
    extract_parametric_coefficients(mediator_models[[med]], exposures) |>
      dplyr::transmute(
        spatial_specification = spatial_spec$name,
        contrast = spatial_contrast_label(term),
        exposure = term,
        mediator = med,
        alpha = estimate,
        SE = std_error
      )
  })

  beta_tbl <- extract_parametric_coefficients(outcome_model, mediators) |>
    dplyr::transmute(
      spatial_specification = spatial_spec$name,
      mediator = term,
      beta = estimate,
      SE = std_error
    )

  direct_tbl <- extract_parametric_coefficients(outcome_model, exposures) |>
    dplyr::transmute(
      spatial_specification = spatial_spec$name,
      contrast = spatial_contrast_label(term),
      exposure = term,
      direct_effect = estimate,
      direct_SE = std_error
    )

  total_tbl <- extract_parametric_coefficients(total_model, exposures) |>
    dplyr::transmute(
      spatial_specification = spatial_spec$name,
      contrast = spatial_contrast_label(term),
      exposure = term,
      TE = estimate,
      TE_SE = std_error
    )

  indirect_tbl <- alpha_tbl |>
    dplyr::left_join(
      beta_tbl |> dplyr::select(spatial_specification, mediator, beta, beta_SE = SE),
      by = c("spatial_specification", "mediator")
    ) |>
    dplyr::mutate(indirect_component = alpha * beta) |>
    dplyr::select(spatial_specification, contrast, exposure, mediator, alpha, beta, indirect_component)

  decomposition_tbl <- indirect_tbl |>
    dplyr::group_by(spatial_specification, contrast, exposure) |>
    dplyr::summarise(total_NIE = sum(indirect_component), .groups = "drop") |>
    dplyr::left_join(total_tbl, by = c("spatial_specification", "contrast", "exposure")) |>
    dplyr::left_join(direct_tbl, by = c("spatial_specification", "contrast", "exposure")) |>
    dplyr::mutate(PM = total_NIE / TE) |>
    dplyr::select(spatial_specification, contrast, exposure, TE, direct_effect, total_NIE, PM, TE_SE, direct_SE)

  list(
    specification = spatial_spec,
    mediator_models = mediator_models,
    outcome_model = outcome_model,
    total_model = total_model,
    residuals = residuals,
    residual_pairwise = make_spatial_pairwise_residual_summary(
      spatial_specification = spatial_spec$name,
      covariance_matrix = residuals$covariance,
      correlation_matrix = residuals$correlation
    ),
    alpha = alpha_tbl,
    beta = beta_tbl,
    direct = direct_tbl,
    total = total_tbl,
    indirect = indirect_tbl,
    decomposition = decomposition_tbl
  )
}

add_change_from_linear <- function(tbl,
                                   value_cols,
                                   by_cols,
                                   spec_col = "spatial_specification",
                                   baseline_spec = "linear") {
  baseline <- tbl |>
    dplyr::filter(.data[[spec_col]] == baseline_spec) |>
    dplyr::select(dplyr::all_of(c(by_cols, value_cols))) |>
    dplyr::rename_with(~ paste0("linear_", .x), dplyr::all_of(value_cols))

  out <- dplyr::left_join(tbl, baseline, by = by_cols)

  for (v in value_cols) {
    linear_v <- paste0("linear_", v)
    out[[paste0(v, "_absolute_change")]] <- out[[v]] - out[[linear_v]]
    out[[paste0(v, "_relative_change")]] <- ifelse(
      abs(out[[linear_v]]) > .Machine$double.eps,
      (out[[v]] - out[[linear_v]]) / out[[linear_v]],
      NA_real_
    )
  }

  out
}

make_knn_neighbor_index <- function(x,
                                    y,
                                    k = 6) {
  coords <- cbind(x, y)
  n <- nrow(coords)
  if (n <= k + 1) {
    stop("Not enough observations for k-nearest-neighbor Moran's I.", call. = FALSE)
  }

  distance <- as.matrix(stats::dist(coords))
  diag(distance) <- Inf
  t(apply(distance, 1, function(row) order(row)[seq_len(k)]))
}

compute_knn_moran_i <- function(value,
                                neighbor_idx) {
  complete <- stats::complete.cases(value)
  value <- value[complete]
  neighbor_idx <- neighbor_idx[complete, , drop = FALSE]
  n <- length(value)
  if (n <= 1) {
    return(NA_real_)
  }

  z <- value - mean(value)
  denom <- sum(z^2)
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }

  numerator <- 0
  for (i in seq_len(n)) {
    numerator <- numerator + sum(z[i] * z[neighbor_idx[i, ]])
  }

  # Row-standardized k-nearest-neighbor weights have total weight n.
  numerator <- numerator / ncol(neighbor_idx)
  (n / n) * numerator / denom
}

compute_spatial_residual_diagnostics <- function(spec_results,
                                                 data,
                                                 mediators = c("PC1_R", "PC2_R", "PC3"),
                                                 outcome = "Y",
                                                 moran_k = 6) {
  neighbor_idx <- make_knn_neighbor_index(data$x_coord, data$y_coord, k = moran_k)

  purrr::imap_dfr(spec_results, function(result, spec_name) {
    mediator_diag <- purrr::map_dfr(mediators, function(med) {
      tibble::tibble(
        spatial_specification = spec_name,
        residual_source = med,
        moran_i_knn = compute_knn_moran_i(
          value = stats::resid(result$mediator_models[[med]], type = "response"),
          neighbor_idx = neighbor_idx
        ),
        knn_k = moran_k
      )
    })

    outcome_diag <- tibble::tibble(
      spatial_specification = spec_name,
      residual_source = outcome,
      moran_i_knn = compute_knn_moran_i(
        value = stats::resid(result$outcome_model, type = "response"),
        neighbor_idx = neighbor_idx
      ),
      knn_k = moran_k
    )

    dplyr::bind_rows(mediator_diag, outcome_diag)
  })
}

run_spatial_sensitivity_analysis <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    coord_x = "imagecol",
    coord_y = "imagerow",
    spline_k = 50,
    spline_basis = "tp",
    include = c("linear", "quadratic", "spline"),
    compute_moran = TRUE,
    moran_k = 6) {
  check_spatial_sensitivity_packages(require_spline = "spline" %in% include)

  analysis_input <- prepare_spatial_sensitivity_data(
    rds_path = rds_path,
    coord_x = coord_x,
    coord_y = coord_y
  )

  specifications <- make_spatial_specifications(
    spline_k = spline_k,
    spline_basis = spline_basis,
    include = include
  )

  spec_results <- purrr::map(
    specifications,
    ~ fit_one_spatial_specification(
      data = analysis_input$data,
      spatial_spec = .x
    )
  )

  residual_correlations <- purrr::map(
    spec_results,
    ~ .x$residuals$correlation
  )
  residual_covariances <- purrr::map(
    spec_results,
    ~ .x$residuals$covariance
  )
  residual_pairwise <- dplyr::bind_rows(purrr::map(spec_results, "residual_pairwise"))

  alpha_comparison <- dplyr::bind_rows(purrr::map(spec_results, "alpha")) |>
    add_change_from_linear(
      value_cols = "alpha",
      by_cols = c("contrast", "exposure", "mediator")
    )

  beta_comparison <- dplyr::bind_rows(purrr::map(spec_results, "beta")) |>
    add_change_from_linear(
      value_cols = "beta",
      by_cols = "mediator"
    )

  indirect_comparison <- dplyr::bind_rows(purrr::map(spec_results, "indirect")) |>
    add_change_from_linear(
      value_cols = c("alpha", "beta", "indirect_component"),
      by_cols = c("contrast", "exposure", "mediator")
    )

  decomposition_comparison <- dplyr::bind_rows(purrr::map(spec_results, "decomposition")) |>
    add_change_from_linear(
      value_cols = c("TE", "direct_effect", "total_NIE", "PM"),
      by_cols = c("contrast", "exposure")
    )

  spatial_diagnostics <- NULL
  if (isTRUE(compute_moran)) {
    spatial_diagnostics <- compute_spatial_residual_diagnostics(
      spec_results = spec_results,
      data = analysis_input$data,
      moran_k = moran_k
    )
  }

  list(
    data = list(
      dat_full_allpc = analysis_input$data,
      pca_df = analysis_input$pca_df,
      pca_fit = analysis_input$pca$pca_fit,
      coordinate_scale_info = analysis_input$coordinate_scale_info,
      predictor_scale_summary = analysis_input$predictor_scale_summary
    ),
    specifications = purrr::map(specifications, function(spec) {
      spec[c("name", "label", "engine", "description", "spline_k", "spline_basis", "method")]
    }),
    model_fits = spec_results,
    residual_covariances = residual_covariances,
    residual_correlations = residual_correlations,
    residual_pairwise = residual_pairwise,
    alpha_comparison = alpha_comparison,
    beta_comparison = beta_comparison,
    indirect_comparison = indirect_comparison,
    decomposition_comparison = decomposition_comparison,
    spatial_diagnostics = spatial_diagnostics,
    notes = c(
      "This is a diagnostic comparison of spatial adjustment specifications, not a final model-selection procedure.",
      "The same spatial specification is applied to mediator, outcome, and total-effect models within each comparison arm.",
      "Quadratic and spline specifications use x_std/y_std built from the current complete-case analysis data.",
      "Because tumor proximity is spatially defined, changes under flexible spatial adjustment should not be interpreted automatically as better confounding control."
    )
  )
}
