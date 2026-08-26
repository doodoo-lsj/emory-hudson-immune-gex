# Frequentist joint mediator diagnostics for the current PCA mediation analysis.
#
# This helper intentionally leaves the main frequentist pipeline unchanged. It
# reuses the current data/PCA preparation and compares separate OLS mediator
# equations with a multivariate Gaussian linear mediator model.

check_joint_mediator_packages <- function() {
  required <- c("dplyr", "tidyr", "tibble", "purrr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Joint mediator diagnostics require installed R packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

contrast_label_from_exposure <- function(exposure) {
  dplyr::case_when(
    exposure == "X_adjacent" ~ "adjacent vs inside",
    exposure == "X_far" ~ "far vs inside",
    TRUE ~ exposure
  )
}

prepare_joint_mediator_data <- function(rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
                                        coord_x = "imagecol",
                                        coord_y = "imagerow",
                                        standardize_covariates = FALSE,
                                        continuous_covariates = c("x_coord", "y_coord")) {
  check_joint_mediator_packages()

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
      ". Source R/pca_mediation_pipeline.R before running joint mediator diagnostics.",
      call. = FALSE
    )
  }

  gex <- load_spatial_mediation_rds(rds_path)
  analysis_df <- prepare_analysis_metadata(gex, coord_x = coord_x, coord_y = coord_y)
  pca <- fit_pca_mediators(gex$M_expr)
  pca_df <- build_pca_score_data(analysis_df, pca$pca_fit)
  dat_full_allpc <- build_full_analysis_data(pca_df)

  dat_fit <- dat_full_allpc
  covariates_for_fit <- continuous_covariates
  scale_info <- NULL

  if (isTRUE(standardize_covariates)) {
    missing_covariates <- setdiff(continuous_covariates, names(dat_full_allpc))
    if (length(missing_covariates) > 0) {
      stop("Missing continuous covariates: ", paste(missing_covariates, collapse = ", "), call. = FALSE)
    }

    scale_rows <- vector("list", length(continuous_covariates))
    for (i in seq_along(continuous_covariates)) {
      v <- continuous_covariates[[i]]
      mu <- mean(dat_full_allpc[[v]], na.rm = TRUE)
      sig <- stats::sd(dat_full_allpc[[v]], na.rm = TRUE)
      if (!is.finite(sig) || sig <= 0) {
        stop("Continuous covariate has zero or invalid SD: ", v, call. = FALSE)
      }

      std_name <- paste0(v, "_std")
      dat_fit[[std_name]] <- (dat_full_allpc[[v]] - mu) / sig

      scale_rows[[i]] <- tibble::tibble(
        variable = v,
        standardized_variable = std_name,
        center = mu,
        scale = sig
      )
    }

    scale_info <- dplyr::bind_rows(scale_rows)
    covariates_for_fit <- scale_info$standardized_variable
  }

  predictor_scale_summary <- purrr::map_dfr(
    c("Y", "PC1_R", "PC2_R", "PC3", continuous_covariates),
    function(v) {
      tibble::tibble(
        variable = v,
        mean = mean(dat_full_allpc[[v]], na.rm = TRUE),
        sd = stats::sd(dat_full_allpc[[v]], na.rm = TRUE),
        min = min(dat_full_allpc[[v]], na.rm = TRUE),
        max = max(dat_full_allpc[[v]], na.rm = TRUE)
      )
    }
  )

  list(
    gex = gex,
    analysis_df = analysis_df,
    pca = pca,
    pca_df = pca_df,
    dat_full_allpc = dat_full_allpc,
    dat_fit = dat_fit,
    covariates_for_fit = covariates_for_fit,
    continuous_covariates = continuous_covariates,
    scale_info = scale_info,
    predictor_scale_summary = predictor_scale_summary
  )
}

build_original_scale_transform <- function(fit_coef_names, scale_info = NULL) {
  original_coef_names <- fit_coef_names

  if (!is.null(scale_info) && nrow(scale_info) > 0) {
    for (i in seq_len(nrow(scale_info))) {
      original_coef_names[original_coef_names == scale_info$standardized_variable[[i]]] <-
        scale_info$variable[[i]]
    }
  }

  transform <- diag(length(fit_coef_names))
  rownames(transform) <- original_coef_names
  colnames(transform) <- fit_coef_names

  if (!is.null(scale_info) && nrow(scale_info) > 0) {
    transform[,] <- 0
    transform["(Intercept)", "(Intercept)"] <- 1

    for (name in fit_coef_names[fit_coef_names != "(Intercept)"]) {
      if (!name %in% scale_info$standardized_variable) {
        transform[name, name] <- 1
      }
    }

    for (i in seq_len(nrow(scale_info))) {
      std_name <- scale_info$standardized_variable[[i]]
      orig_name <- scale_info$variable[[i]]
      transform[orig_name, std_name] <- 1 / scale_info$scale[[i]]
      transform["(Intercept)", std_name] <- -scale_info$center[[i]] / scale_info$scale[[i]]
    }
  }

  transform
}

transform_coefficient_matrix_to_original <- function(coef_matrix, scale_info = NULL) {
  transform <- build_original_scale_transform(rownames(coef_matrix), scale_info)
  transform %*% coef_matrix
}

transform_coefficient_covariance_to_original <- function(coef_covariance,
                                                         fit_coef_names,
                                                         mediator_names,
                                                         scale_info = NULL) {
  transform <- build_original_scale_transform(fit_coef_names, scale_info)
  full_transform <- kronecker(diag(length(mediator_names)), transform)
  original_covariance <- full_transform %*% coef_covariance %*% t(full_transform)

  original_coef_names <- rownames(transform)
  labels <- as.vector(outer(original_coef_names, mediator_names, paste, sep = ":"))
  rownames(original_covariance) <- labels
  colnames(original_covariance) <- labels

  original_covariance
}

tidy_coefficient_matrix <- function(coef_matrix,
                                    se_matrix = NULL,
                                    scale = "original") {
  out <- as.data.frame(coef_matrix) |>
    tibble::rownames_to_column("term") |>
    tidyr::pivot_longer(
      cols = -term,
      names_to = "mediator",
      values_to = "estimate"
    ) |>
    dplyr::mutate(scale = scale)

  if (!is.null(se_matrix)) {
    se_tbl <- as.data.frame(se_matrix) |>
      tibble::rownames_to_column("term") |>
      tidyr::pivot_longer(
        cols = -term,
        names_to = "mediator",
        values_to = "std_error"
      )
    out <- dplyr::left_join(out, se_tbl, by = c("term", "mediator"))
  }

  out
}

fit_separate_mediator_ols <- function(data,
                                      mediators = c("PC1_R", "PC2_R", "PC3"),
                                      exposures = c("X_adjacent", "X_far"),
                                      covariates = c("x_coord", "y_coord"),
                                      scale_info = NULL) {
  predictor_terms <- c(exposures, covariates)
  formulas <- stats::setNames(
    purrr::map(mediators, function(med) {
      stats::as.formula(paste(med, "~", paste(predictor_terms, collapse = " + ")))
    }),
    mediators
  )
  models <- purrr::map(formulas, stats::lm, data = data)

  fit_coef_names <- names(stats::coef(models[[1]]))
  coef_fit <- do.call(
    cbind,
    purrr::map(models, function(fit) stats::coef(fit)[fit_coef_names])
  )
  colnames(coef_fit) <- mediators
  rownames(coef_fit) <- fit_coef_names

  transform <- build_original_scale_transform(fit_coef_names, scale_info)
  coef_original <- transform %*% coef_fit

  se_original <- do.call(
    cbind,
    purrr::map(models, function(fit) {
      vcov_original <- transform %*% stats::vcov(fit)[fit_coef_names, fit_coef_names, drop = FALSE] %*%
        t(transform)
      sqrt(diag(vcov_original))
    })
  )
  colnames(se_original) <- mediators
  rownames(se_original) <- rownames(coef_original)

  residual_matrix <- do.call(cbind, purrr::map(models, stats::resid))
  colnames(residual_matrix) <- mediators
  n <- nrow(residual_matrix)
  p <- length(fit_coef_names)

  sigma_mle <- crossprod(residual_matrix) / n
  sigma_unbiased <- crossprod(residual_matrix) / (n - p)

  list(
    models = models,
    formulas = formulas,
    coefficients_fit_scale = tidy_coefficient_matrix(coef_fit, scale = "fit"),
    coefficients = tidy_coefficient_matrix(coef_original, se_original, scale = "original"),
    residuals = residual_matrix,
    residual_sd = tibble::tibble(
      mediator = mediators,
      residual_sd = sqrt(diag(sigma_unbiased))
    ),
    residual_covariance_mle = sigma_mle,
    residual_covariance_unbiased = sigma_unbiased,
    residual_correlation = stats::cov2cor(sigma_mle)
  )
}

fit_joint_multivariate_mediator <- function(data,
                                            mediators = c("PC1_R", "PC2_R", "PC3"),
                                            exposures = c("X_adjacent", "X_far"),
                                            covariates = c("x_coord", "y_coord"),
                                            scale_info = NULL) {
  predictor_terms <- c(exposures, covariates)
  formula <- stats::as.formula(
    paste0(
      "cbind(",
      paste(mediators, collapse = ", "),
      ") ~ ",
      paste(predictor_terms, collapse = " + ")
    )
  )

  fit <- stats::lm(formula, data = data)
  coef_fit <- stats::coef(fit)
  residual_matrix <- stats::resid(fit)
  colnames(residual_matrix) <- mediators

  x_matrix <- stats::model.matrix(fit)
  n <- nrow(x_matrix)
  p <- ncol(x_matrix)
  xtx_inv <- solve(crossprod(x_matrix))

  sigma_mle <- crossprod(residual_matrix) / n
  sigma_unbiased <- crossprod(residual_matrix) / (n - p)

  coef_cov_fit <- kronecker(sigma_unbiased, xtx_inv)
  fit_labels <- as.vector(outer(rownames(coef_fit), mediators, paste, sep = ":"))
  rownames(coef_cov_fit) <- fit_labels
  colnames(coef_cov_fit) <- fit_labels

  coef_original <- transform_coefficient_matrix_to_original(coef_fit, scale_info)
  coef_cov_original <- transform_coefficient_covariance_to_original(
    coef_covariance = coef_cov_fit,
    fit_coef_names = rownames(coef_fit),
    mediator_names = mediators,
    scale_info = scale_info
  )

  se_original <- matrix(
    sqrt(diag(coef_cov_original)),
    nrow = nrow(coef_original),
    ncol = length(mediators),
    dimnames = list(rownames(coef_original), mediators)
  )

  list(
    fit = fit,
    formula = formula,
    coefficients_fit_scale = tidy_coefficient_matrix(coef_fit, scale = "fit"),
    coefficients = tidy_coefficient_matrix(coef_original, se_original, scale = "original"),
    residuals = residual_matrix,
    residual_covariance_mle = sigma_mle,
    residual_covariance_unbiased = sigma_unbiased,
    residual_correlation = stats::cov2cor(sigma_mle),
    coefficient_covariance_fit_scale = coef_cov_fit,
    coefficient_covariance = coef_cov_original,
    estimator = paste(
      "Coefficient MLE is ordinary least squares for the shared-design",
      "multivariate Gaussian linear model; Sigma_M MLE is E'E / n.",
      "Coefficient covariance is estimated as Sigma_unbiased kron (X'X)^-1."
    )
  )
}

make_residual_pairwise_summary <- function(covariance_matrix, correlation_matrix) {
  mediator_pairs <- utils::combn(colnames(correlation_matrix), 2, simplify = FALSE)

  purrr::map_dfr(mediator_pairs, function(pair) {
    tibble::tibble(
      mediator_1 = pair[[1]],
      mediator_2 = pair[[2]],
      residual_covariance = covariance_matrix[pair[[1]], pair[[2]]],
      residual_correlation = correlation_matrix[pair[[1]], pair[[2]]]
    )
  })
}

make_coefficient_comparison <- function(separate_coefficients,
                                        joint_coefficients,
                                        terms = c("X_adjacent", "X_far", "x_coord", "y_coord")) {
  separate_tbl <- separate_coefficients |>
    dplyr::filter(term %in% terms) |>
    dplyr::select(term, mediator, separate_estimate = estimate, separate_se = std_error)

  joint_tbl <- joint_coefficients |>
    dplyr::filter(term %in% terms) |>
    dplyr::select(term, mediator, joint_estimate = estimate, joint_se = std_error)

  dplyr::left_join(separate_tbl, joint_tbl, by = c("term", "mediator")) |>
    dplyr::mutate(
      difference = joint_estimate - separate_estimate,
      contrast = dplyr::if_else(term %in% c("X_adjacent", "X_far"), contrast_label_from_exposure(term), NA_character_)
    ) |>
    dplyr::select(contrast, mediator, term, separate_estimate, joint_estimate, difference, separate_se, joint_se)
}

extract_alpha_cross_equation_covariance <- function(coefficient_covariance,
                                                    exposures = c("X_adjacent", "X_far"),
                                                    mediators = c("PC1_R", "PC2_R", "PC3")) {
  purrr::map(exposures, function(exposure) {
    labels <- paste(exposure, mediators, sep = ":")
    covariance <- coefficient_covariance[labels, labels, drop = FALSE]
    correlation <- stats::cov2cor(covariance)
    dimnames(covariance) <- list(mediators, mediators)
    dimnames(correlation) <- list(mediators, mediators)

    mediator_pairs <- utils::combn(mediators, 2, simplify = FALSE)
    pairwise <- purrr::map_dfr(mediator_pairs, function(pair) {
      tibble::tibble(
        exposure = exposure,
        contrast = contrast_label_from_exposure(exposure),
        mediator_1 = pair[[1]],
        mediator_2 = pair[[2]],
        alpha_covariance = covariance[pair[[1]], pair[[2]]],
        alpha_correlation = correlation[pair[[1]], pair[[2]]]
      )
    })

    list(
      covariance = covariance,
      correlation = correlation,
      pairwise = pairwise
    )
  }) |>
    stats::setNames(exposures)
}

compute_joint_mediation_from_alpha <- function(data,
                                               joint_coefficients,
                                               exposures = c("X_adjacent", "X_far"),
                                               outcome = "Y",
                                               mediators = c("PC1_R", "PC2_R", "PC3"),
                                               covariates = c("x_coord", "y_coord")) {
  if (!exists("decompose_linear_multix", mode = "function")) {
    stop("decompose_linear_multix() is missing. Source R/pca_mediation_pipeline.R first.", call. = FALSE)
  }

  original <- decompose_linear_multix(
    data = data,
    exposures = exposures,
    outcome = outcome,
    mediators = mediators,
    covariates = covariates
  )

  beta_tbl <- tibble::tibble(
    mediator = mediators,
    beta_M_to_Y = stats::coef(original$fit_outcome)[mediators]
  )

  alpha_tbl <- joint_coefficients |>
    dplyr::filter(term %in% exposures) |>
    dplyr::transmute(
      mediator,
      exposure = term,
      alpha_X_to_M = estimate
    )

  path <- alpha_tbl |>
    dplyr::left_join(beta_tbl, by = "mediator") |>
    dplyr::mutate(
      indirect_component = alpha_X_to_M * beta_M_to_Y,
      contrast = contrast_label_from_exposure(exposure)
    ) |>
    dplyr::select(contrast, exposure, mediator, alpha_X_to_M, beta_M_to_Y, indirect_component)

  summary <- path |>
    dplyr::group_by(exposure, contrast) |>
    dplyr::summarise(NIE = sum(indirect_component), .groups = "drop") |>
    dplyr::mutate(
      TE = stats::coef(original$fit_total)[exposure],
      NDE = stats::coef(original$fit_outcome)[exposure],
      PM = NIE / TE
    ) |>
    dplyr::select(contrast, exposure, TE, NDE, NIE, PM)

  list(
    original = original,
    joint_path = path,
    joint_summary = summary,
    fit_outcome = original$fit_outcome,
    fit_total = original$fit_total
  )
}

make_mediation_comparison <- function(original_decomposition,
                                      joint_path,
                                      joint_summary) {
  original_path <- original_decomposition$path |>
    dplyr::mutate(contrast = contrast_label_from_exposure(exposure)) |>
    dplyr::select(
      contrast,
      exposure,
      mediator,
      original_alpha = alpha_X_to_M,
      beta_M_to_Y,
      original_indirect = indirect_component
    )

  path_comparison <- original_path |>
    dplyr::left_join(
      joint_path |>
        dplyr::select(
          contrast,
          exposure,
          mediator,
          joint_alpha = alpha_X_to_M,
          joint_model_indirect = indirect_component
        ),
      by = c("contrast", "exposure", "mediator")
    ) |>
    dplyr::mutate(
      alpha_difference = joint_alpha - original_alpha,
      indirect_difference = joint_model_indirect - original_indirect
    )

  original_summary <- original_decomposition$summary |>
    dplyr::mutate(contrast = contrast_label_from_exposure(exposure)) |>
    dplyr::select(
      contrast,
      exposure,
      original_TE = TE,
      original_NDE = NDE,
      original_NIE = NIE,
      original_PM = PM
    )

  summary_comparison <- original_summary |>
    dplyr::left_join(
      joint_summary |>
        dplyr::select(
          contrast,
          exposure,
          joint_TE = TE,
          joint_NDE = NDE,
          joint_model_NIE = NIE,
          joint_PM = PM
        ),
      by = c("contrast", "exposure")
    ) |>
    dplyr::mutate(
      TE_difference = joint_TE - original_TE,
      NDE_difference = joint_NDE - original_NDE,
      NIE_difference = joint_model_NIE - original_NIE,
      PM_difference = joint_PM - original_PM
    )

  list(path = path_comparison, summary = summary_comparison)
}

fit_joint_mediator_frequentist_once <- function(rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
                                                coord_x = "imagecol",
                                                coord_y = "imagerow",
                                                standardize_covariates = FALSE,
                                                continuous_covariates = c("x_coord", "y_coord"),
                                                mediators = c("PC1_R", "PC2_R", "PC3"),
                                                exposures = c("X_adjacent", "X_far"),
                                                outcome = "Y",
                                                tolerance = 1e-8) {
  analysis_input <- prepare_joint_mediator_data(
    rds_path = rds_path,
    coord_x = coord_x,
    coord_y = coord_y,
    standardize_covariates = standardize_covariates,
    continuous_covariates = continuous_covariates
  )

  separate <- fit_separate_mediator_ols(
    data = analysis_input$dat_fit,
    mediators = mediators,
    exposures = exposures,
    covariates = analysis_input$covariates_for_fit,
    scale_info = analysis_input$scale_info
  )

  joint <- fit_joint_multivariate_mediator(
    data = analysis_input$dat_fit,
    mediators = mediators,
    exposures = exposures,
    covariates = analysis_input$covariates_for_fit,
    scale_info = analysis_input$scale_info
  )

  coefficient_comparison <- make_coefficient_comparison(
    separate_coefficients = separate$coefficients,
    joint_coefficients = joint$coefficients,
    terms = c(exposures, continuous_covariates)
  )

  residual_pairwise_summary <- make_residual_pairwise_summary(
    covariance_matrix = joint$residual_covariance_mle,
    correlation_matrix = joint$residual_correlation
  )

  mediation <- compute_joint_mediation_from_alpha(
    data = analysis_input$dat_full_allpc,
    joint_coefficients = joint$coefficients,
    exposures = exposures,
    outcome = outcome,
    mediators = mediators,
    covariates = continuous_covariates
  )

  mediation_comparison <- make_mediation_comparison(
    original_decomposition = mediation$original,
    joint_path = mediation$joint_path,
    joint_summary = mediation$joint_summary
  )

  alpha_cross_covariance <- extract_alpha_cross_equation_covariance(
    coefficient_covariance = joint$coefficient_covariance,
    exposures = exposures,
    mediators = mediators
  )

  coefficient_equivalence_max_abs_diff <- max(abs(coefficient_comparison$difference), na.rm = TRUE)
  residual_correlation_max_abs_diff <- max(
    abs(separate$residual_correlation - joint$residual_correlation),
    na.rm = TRUE
  )
  mediation_indirect_max_abs_diff <- max(
    abs(mediation_comparison$path$indirect_difference),
    na.rm = TRUE
  )

  validation <- list(
    coefficient_equivalence_passed = coefficient_equivalence_max_abs_diff <= tolerance,
    coefficient_equivalence_max_abs_diff = coefficient_equivalence_max_abs_diff,
    residual_correlation_equivalence_passed = residual_correlation_max_abs_diff <= tolerance,
    residual_correlation_max_abs_diff = residual_correlation_max_abs_diff,
    mediation_indirect_equivalence_passed = mediation_indirect_max_abs_diff <= tolerance,
    mediation_indirect_max_abs_diff = mediation_indirect_max_abs_diff,
    tolerance = tolerance
  )

  list(
    config = list(
      rds_path = rds_path,
      coord_x = coord_x,
      coord_y = coord_y,
      standardize_covariates = standardize_covariates,
      continuous_covariates = continuous_covariates,
      mediators = mediators,
      exposures = exposures,
      outcome = outcome,
      tolerance = tolerance
    ),
    data = list(
      dat_full_allpc = analysis_input$dat_full_allpc,
      dat_fit = analysis_input$dat_fit,
      predictor_scale_summary = analysis_input$predictor_scale_summary,
      scale_info = analysis_input$scale_info
    ),
    separate_models = separate$models,
    separate = separate,
    joint_model = joint$fit,
    joint = joint,
    coefficient_comparison = coefficient_comparison,
    residual_covariance = joint$residual_covariance_mle,
    residual_covariance_unbiased = joint$residual_covariance_unbiased,
    residual_correlation = joint$residual_correlation,
    residual_pairwise_summary = residual_pairwise_summary,
    coefficient_covariance = joint$coefficient_covariance,
    alpha_cross_equation_covariance = alpha_cross_covariance,
    mediation = mediation,
    mediation_comparison = mediation_comparison,
    validation = validation
  )
}

compare_standardized_joint_results <- function(unstandardized_results,
                                               standardized_results,
                                               tolerance = 1e-8) {
  path_comparison <- unstandardized_results$mediation$joint_path |>
    dplyr::select(contrast, exposure, mediator, unstandardized_indirect = indirect_component) |>
    dplyr::left_join(
      standardized_results$mediation$joint_path |>
        dplyr::select(contrast, exposure, mediator, standardized_indirect = indirect_component),
      by = c("contrast", "exposure", "mediator")
    ) |>
    dplyr::mutate(difference = standardized_indirect - unstandardized_indirect)

  summary_comparison <- unstandardized_results$mediation$joint_summary |>
    dplyr::select(contrast, exposure, unstandardized_NIE = NIE) |>
    dplyr::left_join(
      standardized_results$mediation$joint_summary |>
        dplyr::select(contrast, exposure, standardized_NIE = NIE),
      by = c("contrast", "exposure")
    ) |>
    dplyr::mutate(difference = standardized_NIE - unstandardized_NIE)

  list(
    path = path_comparison,
    summary = summary_comparison,
    passed = max(abs(c(path_comparison$difference, summary_comparison$difference)), na.rm = TRUE) <= tolerance,
    tolerance = tolerance
  )
}

run_joint_mediator_frequentist_analysis <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    coord_x = "imagecol",
    coord_y = "imagerow",
    standardize_covariates = FALSE,
    continuous_covariates = c("x_coord", "y_coord"),
    run_standardization_check = TRUE,
    tolerance = 1e-8) {
  primary <- fit_joint_mediator_frequentist_once(
    rds_path = rds_path,
    coord_x = coord_x,
    coord_y = coord_y,
    standardize_covariates = standardize_covariates,
    continuous_covariates = continuous_covariates,
    tolerance = tolerance
  )

  standardization_check <- NULL
  if (isTRUE(run_standardization_check)) {
    alternate <- fit_joint_mediator_frequentist_once(
      rds_path = rds_path,
      coord_x = coord_x,
      coord_y = coord_y,
      standardize_covariates = !standardize_covariates,
      continuous_covariates = continuous_covariates,
      tolerance = tolerance
    )

    if (isTRUE(standardize_covariates)) {
      standardization_check <- compare_standardized_joint_results(
        unstandardized_results = alternate,
        standardized_results = primary,
        tolerance = tolerance
      )
    } else {
      standardization_check <- compare_standardized_joint_results(
        unstandardized_results = primary,
        standardized_results = alternate,
        tolerance = tolerance
      )
    }
  }

  primary$standardization_check <- standardization_check
  primary
}
