# Spatial overlap diagnostics for the current spline sensitivity result.
#
# This helper keeps the current analysis data, PCA mediators, and spline
# specification fixed. It diagnoses overlap between tumor-proximity exposure
# dummies and the flexible 2D spatial smooth, without selecting a new final
# spatial model.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

check_spatial_overlap_packages <- function() {
  required <- c("dplyr", "tidyr", "tibble", "purrr", "mgcv", "ggplot2")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Spatial overlap diagnostics require installed R packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

overlap_contrast_label <- function(exposure) {
  dplyr::case_when(
    exposure == "X_adjacent" ~ "adjacent vs inside",
    exposure == "X_far" ~ "far vs inside",
    TRUE ~ exposure
  )
}

prepare_spatial_overlap_data <- function(rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
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
      ". Source R/pca_mediation_pipeline.R before running overlap diagnostics.",
      call. = FALSE
    )
  }

  gex <- load_spatial_mediation_rds(rds_path)
  analysis_df <- prepare_analysis_metadata(gex, coord_x = coord_x, coord_y = coord_y)
  pca <- fit_pca_mediators(gex$M_expr)
  pca_df <- build_pca_score_data(analysis_df, pca$pca_fit)
  dat <- build_full_analysis_data(pca_df)

  coord_scale <- tibble::tibble(
    variable = c("x_coord", "y_coord"),
    standardized_variable = c("x_std", "y_std"),
    center = c(mean(dat$x_coord), mean(dat$y_coord)),
    scale = c(stats::sd(dat$x_coord), stats::sd(dat$y_coord))
  )

  if (any(!is.finite(coord_scale$scale)) || any(coord_scale$scale <= 0)) {
    stop("Coordinate SD is zero or invalid.", call. = FALSE)
  }

  dat <- dplyr::mutate(
    dat,
    x_std = (.data$x_coord - coord_scale$center[coord_scale$variable == "x_coord"]) /
      coord_scale$scale[coord_scale$variable == "x_coord"],
    y_std = (.data$y_coord - coord_scale$center[coord_scale$variable == "y_coord"]) /
      coord_scale$scale[coord_scale$variable == "y_coord"]
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
    data = dat,
    gex = gex,
    analysis_df = analysis_df,
    pca = pca,
    pca_df = pca_df,
    coordinate_scale_info = coord_scale,
    predictor_scale_summary = predictor_scale_summary
  )
}

make_overlap_smooth_term <- function(spline_k = 50,
                                     spline_basis = "tp") {
  paste0("s(x_std, y_std, bs = '", spline_basis, "', k = ", spline_k, ")")
}

fit_overlap_gam <- function(response,
                            parametric_terms,
                            data,
                            spline_k = 50,
                            spline_basis = "tp",
                            family = stats::gaussian()) {
  smooth_term <- make_overlap_smooth_term(spline_k = spline_k, spline_basis = spline_basis)
  rhs <- paste(c(parametric_terms, smooth_term), collapse = " + ")
  formula <- stats::as.formula(paste(response, "~", rhs))

  mgcv::gam(formula, data = data, method = "REML", family = family)
}

fit_overlap_gam_set <- function(data,
                                spline_k = 50,
                                spline_basis = "tp",
                                exposures = c("X_adjacent", "X_far"),
                                mediators = c("PC1_R", "PC2_R", "PC3"),
                                outcome = "Y") {
  mediator_models <- stats::setNames(
    purrr::map(mediators, function(med) {
      fit_overlap_gam(
        response = med,
        parametric_terms = exposures,
        data = data,
        spline_k = spline_k,
        spline_basis = spline_basis
      )
    }),
    mediators
  )

  outcome_model <- fit_overlap_gam(
    response = outcome,
    parametric_terms = c(exposures, mediators),
    data = data,
    spline_k = spline_k,
    spline_basis = spline_basis
  )

  total_model <- fit_overlap_gam(
    response = outcome,
    parametric_terms = exposures,
    data = data,
    spline_k = spline_k,
    spline_basis = spline_basis
  )

  list(
    mediator_models = mediator_models,
    outcome_model = outcome_model,
    total_model = total_model
  )
}

fit_linear_mediator_models_for_overlap <- function(data,
                                                   exposures = c("X_adjacent", "X_far"),
                                                   mediators = c("PC1_R", "PC2_R", "PC3")) {
  stats::setNames(
    purrr::map(mediators, function(med) {
      stats::lm(
        stats::as.formula(paste(med, "~", paste(c(exposures, "x_coord", "y_coord"), collapse = " + "))),
        data = data
      )
    }),
    mediators
  )
}

fit_exposure_predictability_gams <- function(data,
                                             spline_k = 50,
                                             spline_basis = "tp",
                                             exposures = c("X_adjacent", "X_far")) {
  stats::setNames(
    purrr::map(exposures, function(exposure) {
      fit_overlap_gam(
        response = exposure,
        parametric_terms = character(0),
        data = data,
        spline_k = spline_k,
        spline_basis = spline_basis,
        family = stats::binomial()
      )
    }),
    exposures
  )
}

extract_gam_fit_summary <- function(fit,
                                    model_name,
                                    model_role,
                                    response) {
  sm <- summary(fit)
  tibble::tibble(
    model = model_name,
    model_role = model_role,
    response = response,
    deviance_explained = unname(sm$dev.expl %||% NA_real_),
    adjusted_r_squared = unname(sm$r.sq %||% NA_real_),
    reml_score = unname(fit$gcv.ubre %||% NA_real_),
    scale_estimate = unname(fit$scale %||% NA_real_),
    n = stats::nobs(fit)
  )
}

extract_smooth_edf <- function(fit,
                               model_name,
                               model_role,
                               response) {
  sm <- summary(fit)
  s_table <- sm$s.table
  if (is.null(s_table) || nrow(s_table) == 0) {
    return(tibble::tibble(
      model = model_name,
      model_role = model_role,
      response = response,
      smooth = NA_character_,
      edf = NA_real_,
      reference_df = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_
    ))
  }

  statistic_col <- intersect(c("F", "Chi.sq"), colnames(s_table))[1] %||% NA_character_
  p_col <- intersect(c("p-value", "p.value"), colnames(s_table))[1] %||% NA_character_

  tibble::tibble(
    model = model_name,
    model_role = model_role,
    response = response,
    smooth = rownames(s_table),
    edf = s_table[, "edf"],
    reference_df = s_table[, "Ref.df"],
    statistic = if (!is.na(statistic_col)) s_table[, statistic_col] else NA_real_,
    p_value = if (!is.na(p_col)) s_table[, p_col] else NA_real_
  )
}

tidy_concurvity_full <- function(concurvity_object,
                                 model_name) {
  if (!is.matrix(concurvity_object) && !is.data.frame(concurvity_object)) {
    return(tibble::tibble())
  }

  as.data.frame(concurvity_object) |>
    tibble::rownames_to_column("metric") |>
    tidyr::pivot_longer(
      cols = -metric,
      names_to = "term",
      values_to = "concurvity"
    ) |>
    dplyr::mutate(model = model_name, .before = 1)
}

tidy_concurvity_pairwise <- function(concurvity_object,
                                     model_name) {
  if (!is.list(concurvity_object)) {
    return(tibble::tibble())
  }

  purrr::imap_dfr(concurvity_object, function(mat, metric) {
    if (!is.matrix(mat) && !is.data.frame(mat)) {
      return(tibble::tibble())
    }

    as.data.frame(mat) |>
      tibble::rownames_to_column("term_from") |>
      tidyr::pivot_longer(
        cols = -term_from,
        names_to = "term_to",
        values_to = "concurvity"
      ) |>
      dplyr::mutate(model = model_name, metric = metric, .before = 1)
  })
}

compute_concurvity_diagnostics <- function(gam_fits) {
  all_models <- c(
    gam_fits$mediator_models,
    list(outcome = gam_fits$outcome_model, total = gam_fits$total_model)
  )

  raw_full <- purrr::map(all_models, ~ mgcv::concurvity(.x, full = TRUE))
  raw_pairwise <- purrr::map(all_models, ~ mgcv::concurvity(.x, full = FALSE))

  list(
    raw_full = raw_full,
    raw_pairwise = raw_pairwise,
    tidy_full = purrr::imap_dfr(raw_full, tidy_concurvity_full),
    tidy_pairwise = purrr::imap_dfr(raw_pairwise, tidy_concurvity_pairwise),
    note = paste(
      "mgcv concurvity is returned in both raw and tidy forms.",
      "If parametric dummy terms are grouped as 'para', use the exposure spatial",
      "predictability diagnostics as the direct X-vs-coordinate overlap measure."
    )
  )
}

compute_auc <- function(y,
                        score) {
  y <- as.integer(y)
  complete <- stats::complete.cases(y, score)
  y <- y[complete]
  score <- score[complete]

  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  if (n_pos == 0 || n_neg == 0) {
    return(NA_real_)
  }

  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

summarize_probability_distribution <- function(data,
                                               exposure,
                                               predicted_probability) {
  tibble::tibble(
    X_cat = data$X_cat,
    observed = data[[exposure]],
    predicted_probability = predicted_probability
  ) |>
    dplyr::group_by(X_cat, observed) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(predicted_probability),
      sd = stats::sd(predicted_probability),
      min = min(predicted_probability),
      q025 = stats::quantile(predicted_probability, 0.025),
      median = stats::median(predicted_probability),
      q975 = stats::quantile(predicted_probability, 0.975),
      max = max(predicted_probability),
      .groups = "drop"
    ) |>
    dplyr::mutate(exposure = exposure, .before = 1)
}

compute_exposure_spatial_predictability <- function(exposure_models,
                                                    data,
                                                    exposures = c("X_adjacent", "X_far")) {
  summaries <- purrr::map_dfr(exposures, function(exposure) {
    fit <- exposure_models[[exposure]]
    predicted <- as.numeric(stats::predict(fit, type = "response"))
    sm <- summary(fit)
    observed <- data[[exposure]]

    tibble::tibble(
      exposure = exposure,
      contrast = overlap_contrast_label(exposure),
      auc = compute_auc(observed, predicted),
      accuracy_at_0_5 = mean(as.integer(predicted >= 0.5) == observed),
      deviance_explained = unname(sm$dev.expl %||% NA_real_),
      adjusted_r_squared = unname(sm$r.sq %||% NA_real_),
      mean_predicted_probability = mean(predicted),
      sd_predicted_probability = stats::sd(predicted),
      min_predicted_probability = min(predicted),
      median_predicted_probability = stats::median(predicted),
      max_predicted_probability = max(predicted),
      event_rate = mean(observed)
    )
  })

  probability_distribution <- purrr::map_dfr(exposures, function(exposure) {
    summarize_probability_distribution(
      data = data,
      exposure = exposure,
      predicted_probability = as.numeric(stats::predict(exposure_models[[exposure]], type = "response"))
    )
  })

  list(
    fits = exposure_models,
    summary = summaries,
    probability_distribution = probability_distribution
  )
}

get_spatial_smooth_term_column <- function(term_matrix) {
  term_cols <- colnames(term_matrix)
  smooth_cols <- grep("^s\\(x_std,y_std\\)|^s\\(x_std, y_std\\)", term_cols, value = TRUE)
  if (length(smooth_cols) == 0) {
    smooth_cols <- grep("x_std.*y_std|y_std.*x_std", term_cols, value = TRUE)
  }
  if (length(smooth_cols) == 0) {
    stop("Could not identify spatial smooth contribution column from predict(type = 'terms').", call. = FALSE)
  }
  smooth_cols[[1]]
}

summarize_smooth_by_exposure <- function(gam_fits,
                                         data,
                                         mediators = c("PC1_R", "PC2_R", "PC3")) {
  smooth_summary <- purrr::map_dfr(mediators, function(med) {
    term_matrix <- stats::predict(gam_fits$mediator_models[[med]], type = "terms")
    smooth_col <- get_spatial_smooth_term_column(term_matrix)

    tibble::tibble(
      mediator = med,
      X_cat = data$X_cat,
      smooth_contribution = as.numeric(term_matrix[, smooth_col])
    ) |>
      dplyr::group_by(mediator, X_cat) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean = mean(smooth_contribution),
        sd = stats::sd(smooth_contribution),
        min = min(smooth_contribution),
        median = stats::median(smooth_contribution),
        max = max(smooth_contribution),
        .groups = "drop"
      )
  })

  smooth_wide <- smooth_summary |>
    dplyr::select(mediator, X_cat, mean) |>
    tidyr::pivot_wider(names_from = X_cat, values_from = mean)

  smooth_differences <- smooth_wide |>
    dplyr::transmute(
      mediator,
      adjacent_minus_inside = adjacent - inside,
      far_minus_inside = far - inside,
      far_minus_adjacent = far - adjacent
    )

  list(
    group_summary = smooth_summary,
    group_mean_differences = smooth_differences
  )
}

extract_alpha_terms <- function(model,
                                mediator,
                                model_type,
                                exposures = c("X_adjacent", "X_far")) {
  coef_table <- if (inherits(model, "gam")) {
    summary(model)$p.table
  } else {
    summary(model)$coefficients
  }

  tibble::tibble(
    mediator = mediator,
    model_type = model_type,
    exposure = exposures,
    contrast = overlap_contrast_label(exposures),
    estimate = coef_table[exposures, "Estimate"],
    SE = coef_table[exposures, "Std. Error"],
    abs_estimate_over_SE = abs(estimate) / SE
  )
}

compute_alpha_stability <- function(linear_mediator_models,
                                    spline_mediator_models,
                                    exposures = c("X_adjacent", "X_far"),
                                    mediators = c("PC1_R", "PC2_R", "PC3")) {
  linear_alpha <- purrr::map_dfr(mediators, function(med) {
    extract_alpha_terms(linear_mediator_models[[med]], med, "linear", exposures)
  })
  spline_alpha <- purrr::map_dfr(mediators, function(med) {
    extract_alpha_terms(spline_mediator_models[[med]], med, "spline", exposures)
  })

  linear_alpha |>
    dplyr::select(
      mediator,
      exposure,
      contrast,
      linear_alpha = estimate,
      linear_SE = SE,
      linear_abs_estimate_over_SE = abs_estimate_over_SE
    ) |>
    dplyr::left_join(
      spline_alpha |>
        dplyr::select(
          mediator,
          exposure,
          spline_alpha = estimate,
          spline_SE = SE,
          spline_abs_estimate_over_SE = abs_estimate_over_SE
        ),
      by = c("mediator", "exposure")
    ) |>
    dplyr::mutate(
      absolute_change = spline_alpha - linear_alpha,
      percent_change = ifelse(
        abs(linear_alpha) > .Machine$double.eps,
        100 * (spline_alpha - linear_alpha) / linear_alpha,
        NA_real_
      ),
      SE_ratio = spline_SE / linear_SE
    ) |>
    dplyr::select(
      contrast,
      mediator,
      exposure,
      linear_alpha,
      spline_alpha,
      absolute_change,
      percent_change,
      linear_SE,
      spline_SE,
      SE_ratio,
      linear_abs_estimate_over_SE,
      spline_abs_estimate_over_SE
    )
}

extract_k_diagnostics <- function(fit,
                                  model_name,
                                  model_role,
                                  response) {
  k_check_fun <- getFromNamespace("k.check", "mgcv")
  k_result <- tryCatch(k_check_fun(fit), error = function(e) e)

  if (inherits(k_result, "error") || is.null(k_result)) {
    return(tibble::tibble(
      model = model_name,
      model_role = model_role,
      response = response,
      smooth = NA_character_,
      k_prime = NA_real_,
      edf = NA_real_,
      k_index = NA_real_,
      p_value = NA_real_,
      note = paste("k.check failed:", conditionMessage(k_result))
    ))
  }

  k_df <- as.data.frame(k_result) |>
    tibble::rownames_to_column("smooth")

  names(k_df) <- gsub("k'", "k_prime", names(k_df), fixed = TRUE)
  names(k_df) <- gsub("k-index", "k_index", names(k_df), fixed = TRUE)
  names(k_df) <- gsub("p-value", "p_value", names(k_df), fixed = TRUE)

  k_df |>
    dplyr::mutate(
      model = model_name,
      model_role = model_role,
      response = response,
      note = NA_character_,
      .before = 1
    )
}

make_exposure_plots <- function(data,
                                exposure_predictability,
                                coordinate_scale_info,
                                grid_n = 120) {
  exposure_map <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = x_coord, y = y_coord, color = X_cat)
  ) +
    ggplot2::geom_point(size = 1.8, alpha = 0.85) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "x coordinate", y = "y coordinate", color = "X category")

  x_grid <- seq(min(data$x_std), max(data$x_std), length.out = grid_n)
  y_grid <- seq(min(data$y_std), max(data$y_std), length.out = grid_n)
  pred_grid <- expand.grid(x_std = x_grid, y_std = y_grid)

  x_info <- coordinate_scale_info[coordinate_scale_info$variable == "x_coord", ]
  y_info <- coordinate_scale_info[coordinate_scale_info$variable == "y_coord", ]
  pred_grid$x_coord <- pred_grid$x_std * x_info$scale + x_info$center
  pred_grid$y_coord <- pred_grid$y_std * y_info$scale + y_info$center

  make_probability_plot <- function(exposure, title) {
    fit <- exposure_predictability$fits[[exposure]]
    plot_df <- pred_grid
    plot_df$predicted_probability <- as.numeric(stats::predict(fit, newdata = pred_grid, type = "response"))

    ggplot2::ggplot(plot_df, ggplot2::aes(x = x_coord, y = y_coord, fill = predicted_probability)) +
      ggplot2::geom_raster() +
      ggplot2::geom_point(
        data = data,
        ggplot2::aes(x = x_coord, y = y_coord),
        inherit.aes = FALSE,
        size = 0.25,
        alpha = 0.25
      ) +
      ggplot2::coord_equal() +
      ggplot2::scale_fill_viridis_c(limits = c(0, 1)) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "x coordinate", y = "y coordinate", fill = "Pr", title = title)
  }

  list(
    exposure_map = exposure_map,
    adjacent_probability_surface = make_probability_plot(
      "X_adjacent",
      "P(X_adjacent = 1 | x, y)"
    ),
    far_probability_surface = make_probability_plot(
      "X_far",
      "P(X_far = 1 | x, y)"
    )
  )
}

run_spatial_overlap_diagnostics <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    coord_x = "imagecol",
    coord_y = "imagerow",
    spline_k = 50,
    spline_basis = "tp",
    grid_n = 120) {
  check_spatial_overlap_packages()

  analysis_input <- prepare_spatial_overlap_data(
    rds_path = rds_path,
    coord_x = coord_x,
    coord_y = coord_y
  )
  data <- analysis_input$data

  exposures <- c("X_adjacent", "X_far")
  mediators <- c("PC1_R", "PC2_R", "PC3")

  spline_fits <- fit_overlap_gam_set(
    data = data,
    spline_k = spline_k,
    spline_basis = spline_basis,
    exposures = exposures,
    mediators = mediators,
    outcome = "Y"
  )

  linear_mediator_models <- fit_linear_mediator_models_for_overlap(
    data = data,
    exposures = exposures,
    mediators = mediators
  )

  exposure_models <- fit_exposure_predictability_gams(
    data = data,
    spline_k = spline_k,
    spline_basis = spline_basis,
    exposures = exposures
  )

  fit_index <- dplyr::bind_rows(
    tibble::tibble(
      model = mediators,
      model_role = "mediator",
      response = mediators,
      fit = unname(spline_fits$mediator_models)
    ),
    tibble::tibble(
      model = c("outcome", "total"),
      model_role = c("outcome", "total"),
      response = c("Y", "Y"),
      fit = list(spline_fits$outcome_model, spline_fits$total_model)
    )
  )

  gam_fit_summary <- purrr::pmap_dfr(
    fit_index,
    function(model, model_role, response, fit) {
      extract_gam_fit_summary(fit, model, model_role, response)
    }
  )

  smooth_edf <- purrr::pmap_dfr(
    fit_index,
    function(model, model_role, response, fit) {
      extract_smooth_edf(fit, model, model_role, response)
    }
  )

  k_diagnostics <- purrr::pmap_dfr(
    fit_index,
    function(model, model_role, response, fit) {
      extract_k_diagnostics(fit, model, model_role, response)
    }
  )

  concurvity <- compute_concurvity_diagnostics(spline_fits)
  exposure_spatial_predictability <- compute_exposure_spatial_predictability(
    exposure_models = exposure_models,
    data = data,
    exposures = exposures
  )
  smooth_by_exposure <- summarize_smooth_by_exposure(
    gam_fits = spline_fits,
    data = data,
    mediators = mediators
  )
  alpha_stability <- compute_alpha_stability(
    linear_mediator_models = linear_mediator_models,
    spline_mediator_models = spline_fits$mediator_models,
    exposures = exposures,
    mediators = mediators
  )
  plots <- make_exposure_plots(
    data = data,
    exposure_predictability = exposure_spatial_predictability,
    coordinate_scale_info = analysis_input$coordinate_scale_info,
    grid_n = grid_n
  )

  list(
    config = list(
      rds_path = rds_path,
      coord_x = coord_x,
      coord_y = coord_y,
      spline_k = spline_k,
      spline_basis = spline_basis,
      spline_specification = make_overlap_smooth_term(spline_k, spline_basis),
      method = "REML"
    ),
    data = list(
      dat_full_allpc = data,
      pca_df = analysis_input$pca_df,
      pca_fit = analysis_input$pca$pca_fit,
      coordinate_scale_info = analysis_input$coordinate_scale_info,
      predictor_scale_summary = analysis_input$predictor_scale_summary
    ),
    fits = list(
      spline = spline_fits,
      linear_mediator_models = linear_mediator_models,
      exposure_spatial_predictability = exposure_models
    ),
    gam_fit_summary = gam_fit_summary,
    smooth_edf = smooth_edf,
    concurvity = concurvity,
    exposure_spatial_predictability = exposure_spatial_predictability,
    smooth_by_exposure = smooth_by_exposure,
    alpha_stability = alpha_stability,
    k_diagnostics = k_diagnostics,
    plots = plots,
    notes = c(
      "This diagnostic quantifies overlap between spatially defined exposure and flexible spatial adjustment.",
      "It does not determine whether the spline model is correct or whether the linear model is correct.",
      "High coordinate-based exposure predictability may reflect the definition of tumor proximity, not necessarily a modeling error.",
      "The spline specification is fixed to mgcv::s(x_std, y_std, bs = 'tp', k = 50), method = 'REML'."
    )
  )
}
