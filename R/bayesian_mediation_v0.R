# Bayesian v0 helper functions for the current PCA mediation analysis.
#
# This file intentionally mirrors the current frequentist linear mediation
# specification. It treats PC1_R, PC2_R, and PC3 as fixed observed mediators
# constructed by R/pca_mediation_pipeline.R. It does not implement Bayesian PCA,
# multivariate mediator models, residual covariance, spatial random effects, or
# sensitivity analysis.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

check_bayesian_v0_packages <- function() {
  required <- c("brms", "posterior", "dplyr", "tidyr", "tibble", "purrr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Bayesian v0 requires installed R packages: ",
      paste(missing, collapse = ", "),
      ". Install them before running scripts/current/run_bayesian_v0.R.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

summarize_bayesian_v0_predictor_scales <- function(data,
                                                   variables = c(
                                                     "Y",
                                                     "PC1_R",
                                                     "PC2_R",
                                                     "PC3",
                                                     "x_coord",
                                                     "y_coord"
                                                   )) {
  missing <- setdiff(variables, names(data))
  if (length(missing) > 0) {
    stop("Missing variables for scale summary: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  purrr::map_dfr(variables, function(v) {
    x <- data[[v]]
    tibble::tibble(
      variable = v,
      mean = mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  })
}

standardize_bayesian_v0_coordinates <- function(data,
                                                coord_vars = c("x_coord", "y_coord")) {
  missing <- setdiff(coord_vars, names(data))
  if (length(missing) > 0) {
    stop("Missing coordinate variables: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  data_std <- data
  scale_info <- vector("list", length(coord_vars))

  for (i in seq_along(coord_vars)) {
    v <- coord_vars[[i]]
    mu <- mean(data[[v]], na.rm = TRUE)
    sig <- stats::sd(data[[v]], na.rm = TRUE)
    if (!is.finite(sig) || sig <= 0) {
      stop("Coordinate has zero or invalid SD: ", v, call. = FALSE)
    }

    data_std[[paste0(v, "_std")]] <- (data[[v]] - mu) / sig

    scale_info[[i]] <- tibble::tibble(
      variable = v,
      standardized_variable = paste0(v, "_std"),
      center = mu,
      scale = sig
    )
  }

  list(data = data_std, coordinate_scale_info = dplyr::bind_rows(scale_info))
}

prepare_bayesian_mediation_v0_data <- function(rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
                                               coord_x = "imagecol",
                                               coord_y = "imagerow",
                                               standardize_coordinates = TRUE) {
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
      ". Source R/pca_mediation_pipeline.R before using Bayesian v0.",
      call. = FALSE
    )
  }

  gex <- load_spatial_mediation_rds(rds_path)
  analysis_df <- prepare_analysis_metadata(gex, coord_x = coord_x, coord_y = coord_y)
  pca <- fit_pca_mediators(gex$M_expr)
  pca_df <- build_pca_score_data(analysis_df, pca$pca_fit)
  dat_full_allpc <- build_full_analysis_data(pca_df)
  predictor_scale_summary <- summarize_bayesian_v0_predictor_scales(dat_full_allpc)

  coordinate_scale_info <- NULL
  dat_bayes <- dat_full_allpc
  if (isTRUE(standardize_coordinates)) {
    standardized <- standardize_bayesian_v0_coordinates(dat_full_allpc)
    dat_bayes <- standardized$data
    coordinate_scale_info <- standardized$coordinate_scale_info
  }

  list(
    gex = gex,
    analysis_df = analysis_df,
    pca = pca,
    pca_df = pca_df,
    dat_full_allpc = dat_full_allpc,
    dat_bayes = dat_bayes,
    predictor_scale_summary = predictor_scale_summary,
    coordinate_scale_info = coordinate_scale_info
  )
}

bayesian_mediation_v0_formulas <- function(use_standardized_coordinates = TRUE) {
  coord_terms <- if (isTRUE(use_standardized_coordinates)) {
    c("x_coord_std", "y_coord_std")
  } else {
    c("x_coord", "y_coord")
  }
  coord_part <- paste(coord_terms, collapse = " + ")

  # These formulas are algebraically equivalent to the current frequentist
  # formulas when standardized coordinates are used. Only the numerical
  # parameterization of the coordinate coefficients changes; exposure and PC
  # coefficients used for mediation remain on the original analysis scale.
  list(
    mediator_PC1_R = stats::as.formula(paste("PC1_R ~ X_adjacent + X_far +", coord_part)),
    mediator_PC2_R = stats::as.formula(paste("PC2_R ~ X_adjacent + X_far +", coord_part)),
    mediator_PC3 = stats::as.formula(paste("PC3 ~ X_adjacent + X_far +", coord_part)),
    outcome = stats::as.formula(paste("Y ~ X_adjacent + X_far + PC1_R + PC2_R + PC3 +", coord_part)),
    total = stats::as.formula(paste("Y ~ X_adjacent + X_far +", coord_part))
  )
}

format_prior_number <- function(x) {
  formatC(x, digits = 8, format = "fg", flag = "#")
}

make_weak_scale_priors_for_formula <- function(data, formula) {
  check_bayesian_v0_packages()

  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  y_sd <- stats::sd(y, na.rm = TRUE)
  y_mean <- mean(y, na.rm = TRUE)

  if (!is.finite(y_sd) || y_sd <= 0) {
    stop("Response has zero or invalid standard deviation for formula: ", deparse(formula), call. = FALSE)
  }

  x <- stats::model.matrix(formula, data = data)
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]

  # Priors are data-scale weak priors. A one-SD predictor change is allowed to
  # move the response by roughly 10 response SDs a priori, avoiding the strong
  # shrinkage that would result from using the same Normal(0, 1) prior for all
  # coefficients when coordinates and PC scores have very different scales.
  priors <- c(
    brms::prior_string(
      paste0(
        "student_t(3, ",
        format_prior_number(y_mean),
        ", ",
        format_prior_number(10 * y_sd),
        ")"
      ),
      class = "Intercept"
    ),
    brms::prior_string(
      paste0("exponential(", format_prior_number(1 / y_sd), ")"),
      class = "sigma"
    )
  )

  for (coef_name in colnames(x)) {
    x_sd <- stats::sd(x[, coef_name], na.rm = TRUE)
    if (!is.finite(x_sd) || x_sd <= 0) {
      stop("Predictor has zero or invalid SD: ", coef_name, " in formula: ", deparse(formula), call. = FALSE)
    }

    coef_scale <- 10 * y_sd / x_sd
    priors <- c(
      priors,
      brms::prior_string(
        paste0("student_t(3, 0, ", format_prior_number(coef_scale), ")"),
        class = "b",
        coef = coef_name
      )
    )
  }

  priors
}

make_bayesian_mediation_v0_priors <- function(data,
                                              formulas = bayesian_mediation_v0_formulas()) {
  purrr::map(formulas, ~ make_weak_scale_priors_for_formula(data, .x))
}

fit_single_bayesian_lm <- function(formula,
                                   data,
                                   prior,
                                   chains = 4,
                                   iter = 2000,
                                   warmup = 1000,
                                   cores = 4,
                                   seed = 123,
                                   backend = NULL,
                                   adapt_delta = NULL,
                                   max_treedepth = NULL,
                                   refresh = 100) {
  check_bayesian_v0_packages()

  control <- list()
  if (!is.null(adapt_delta)) {
    control$adapt_delta <- adapt_delta
  }
  if (!is.null(max_treedepth)) {
    control$max_treedepth <- max_treedepth
  }

  brm_args <- list(
    formula = brms::bf(formula),
    data = data,
    family = stats::gaussian(),
    prior = prior,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    refresh = refresh,
    control = control
  )

  if (!is.null(backend)) {
    brm_args$backend <- backend
  }

  do.call(brms::brm, brm_args)
}

fit_bayesian_mediation_v0 <- function(data,
                                      formulas = bayesian_mediation_v0_formulas(),
                                      priors = make_bayesian_mediation_v0_priors(data, formulas),
                                      chains = 4,
                                      iter = 2000,
                                      warmup = 1000,
                                      cores = 4,
                                      seed = 123,
                                      backend = NULL,
                                      adapt_delta = NULL,
                                      max_treedepth = NULL,
                                      refresh = 100) {
  check_bayesian_v0_packages()

  timing <- list()

  fit_one <- function(model_name, model_seed) {
    elapsed <- system.time({
      fit <- fit_single_bayesian_lm(
        formula = formulas[[model_name]],
        data = data,
        prior = priors[[model_name]],
        chains = chains,
        iter = iter,
        warmup = warmup,
        cores = cores,
        seed = model_seed,
        backend = backend,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        refresh = refresh
      )
    })
    timing[[model_name]] <<- tibble::tibble(
      component = model_name,
      elapsed_seconds = unname(elapsed[["elapsed"]])
    )
    fit
  }

  fits <- list(
    mediator_PC1_R = fit_one("mediator_PC1_R", seed + 1),
    mediator_PC2_R = fit_one("mediator_PC2_R", seed + 2),
    mediator_PC3 = fit_one("mediator_PC3", seed + 3),
    outcome = fit_one("outcome", seed + 4),
    total = fit_one("total", seed + 5)
  )

  attr(fits, "timing") <- dplyr::bind_rows(timing)
  fits
}

extract_b_parameter <- function(draws_df, coef_name) {
  var_name <- paste0("b_", coef_name)
  if (!var_name %in% names(draws_df)) {
    stop("Coefficient not found in posterior draws: ", var_name, call. = FALSE)
  }
  draws_df[[var_name]]
}

posterior_summary_vector <- function(x) {
  tibble::tibble(
    mean = mean(x, na.rm = TRUE),
    median = stats::median(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    q025 = stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    q975 = stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE),
    Pr_gt_0 = mean(x > 0, na.rm = TRUE),
    Pr_lt_0 = mean(x < 0, na.rm = TRUE)
  )
}

summarize_original_scale_fixed_effects <- function(fits,
                                                   coordinate_scale_info = NULL) {
  check_bayesian_v0_packages()

  purrr::imap_dfr(fits, function(fit, model) {
    draws <- posterior::as_draws_df(fit)
    coef_names <- names(draws)[grepl("^b_", names(draws))]

    if (length(coef_names) == 0) {
      return(tibble::tibble())
    }

    coef_draws <- draws[, coef_names, drop = FALSE]
    names(coef_draws) <- sub("^b_", "", names(coef_draws))

    if (!is.null(coordinate_scale_info) && nrow(coordinate_scale_info) > 0) {
      for (i in seq_len(nrow(coordinate_scale_info))) {
        original_name <- coordinate_scale_info$variable[[i]]
        std_name <- coordinate_scale_info$standardized_variable[[i]]
        center <- coordinate_scale_info$center[[i]]
        scale <- coordinate_scale_info$scale[[i]]

        if (std_name %in% names(coef_draws)) {
          if ("Intercept" %in% names(coef_draws)) {
            coef_draws$Intercept <- coef_draws$Intercept - (center / scale) * coef_draws[[std_name]]
          }
          coef_draws[[original_name]] <- coef_draws[[std_name]] / scale
          coef_draws[[std_name]] <- NULL
        }
      }
    }

    purrr::map_dfr(names(coef_draws), function(term) {
      s <- posterior_summary_vector(coef_draws[[term]])
      tibble::tibble(
        model = model,
        term = term,
        mean = s$mean,
        median = s$median,
        sd = s$sd,
        q025 = s$q025,
        q975 = s$q975,
        Pr_gt_0 = s$Pr_gt_0,
        Pr_lt_0 = s$Pr_lt_0
      )
    })
  })
}

summarize_bayesian_paths <- function(posterior_quantities) {
  d <- posterior_quantities$draws

  path_specs <- tibble::tribble(
    ~contrast, ~mediator, ~alpha_col, ~beta_col, ~indirect_col,
    "adjacent vs inside", "PC1_R", "alpha_A_PC1_R", "beta_PC1_R", "IE_PC1_A",
    "adjacent vs inside", "PC2_R", "alpha_A_PC2_R", "beta_PC2_R", "IE_PC2_A",
    "adjacent vs inside", "PC3", "alpha_A_PC3", "beta_PC3", "IE_PC3_A",
    "far vs inside", "PC1_R", "alpha_F_PC1_R", "beta_PC1_R", "IE_PC1_F",
    "far vs inside", "PC2_R", "alpha_F_PC2_R", "beta_PC2_R", "IE_PC2_F",
    "far vs inside", "PC3", "alpha_F_PC3", "beta_PC3", "IE_PC3_F"
  )

  purrr::pmap_dfr(path_specs, function(contrast, mediator, alpha_col, beta_col, indirect_col) {
    alpha_sum <- posterior_summary_vector(d[[alpha_col]])
    beta_sum <- posterior_summary_vector(d[[beta_col]])
    indirect_sum <- posterior_summary_vector(d[[indirect_col]])

    tibble::tibble(
      contrast = contrast,
      mediator = mediator,
      alpha_mean = alpha_sum$mean,
      alpha_q025 = alpha_sum$q025,
      alpha_q975 = alpha_sum$q975,
      beta_mean = beta_sum$mean,
      beta_q025 = beta_sum$q025,
      beta_q975 = beta_sum$q975,
      indirect_mean = indirect_sum$mean,
      indirect_median = indirect_sum$median,
      indirect_q025 = indirect_sum$q025,
      indirect_q975 = indirect_sum$q975,
      Pr_indirect_gt_0 = indirect_sum$Pr_gt_0
    )
  })
}

summarize_bayesian_decomposition <- function(posterior_quantities) {
  d <- posterior_quantities$draws

  decomp_specs <- tibble::tribble(
    ~contrast, ~te_col, ~direct_col, ~nie_col, ~pm_col,
    "adjacent vs inside", "TE_A", "DE_A", "NIE_A", "PM_A",
    "far vs inside", "TE_F", "DE_F", "NIE_F", "PM_F"
  )

  purrr::pmap_dfr(decomp_specs, function(contrast, te_col, direct_col, nie_col, pm_col) {
    te_sum <- posterior_summary_vector(d[[te_col]])
    direct_sum <- posterior_summary_vector(d[[direct_col]])
    nie_sum <- posterior_summary_vector(d[[nie_col]])
    pm_sum <- posterior_summary_vector(d[[pm_col]])

    tibble::tibble(
      contrast = contrast,
      TE_mean = te_sum$mean,
      TE_q025 = te_sum$q025,
      TE_q975 = te_sum$q975,
      direct_mean = direct_sum$mean,
      direct_q025 = direct_sum$q025,
      direct_q975 = direct_sum$q975,
      NIE_mean = nie_sum$mean,
      NIE_median = nie_sum$median,
      NIE_q025 = nie_sum$q025,
      NIE_q975 = nie_sum$q975,
      PM_mean = pm_sum$mean,
      PM_q025 = pm_sum$q025,
      PM_q975 = pm_sum$q975,
      Pr_NIE_gt_0 = nie_sum$Pr_gt_0
    )
  })
}

compute_bayesian_mediation_quantities <- function(fits,
                                                  n_draws = NULL,
                                                  seed = 123) {
  check_bayesian_v0_packages()
  set.seed(seed)

  draws <- list(
    mediator_PC1_R = posterior::as_draws_df(fits$mediator_PC1_R),
    mediator_PC2_R = posterior::as_draws_df(fits$mediator_PC2_R),
    mediator_PC3 = posterior::as_draws_df(fits$mediator_PC3),
    outcome = posterior::as_draws_df(fits$outcome),
    total = posterior::as_draws_df(fits$total)
  )

  available_draws <- vapply(draws, nrow, integer(1))
  if (is.null(n_draws)) {
    n_draws <- min(available_draws)
  }
  if (n_draws <= 0) {
    stop("n_draws must be positive.", call. = FALSE)
  }

  # The brms models are fit independently. Therefore, identical row or
  # iteration numbers across models are not treated as a joint posterior draw.
  # Instead, each posterior Monte Carlo row below independently resamples one
  # draw from each fitted model, preserving within-model coefficient dependence
  # while representing the product posterior under independent model fits.
  idx <- lapply(available_draws, function(n) sample.int(n, size = n_draws, replace = TRUE))
  names(idx) <- names(available_draws)

  d_pc1 <- draws$mediator_PC1_R[idx$mediator_PC1_R, , drop = FALSE]
  d_pc2 <- draws$mediator_PC2_R[idx$mediator_PC2_R, , drop = FALSE]
  d_pc3 <- draws$mediator_PC3[idx$mediator_PC3, , drop = FALSE]
  d_out <- draws$outcome[idx$outcome, , drop = FALSE]
  d_total <- draws$total[idx$total, , drop = FALSE]

  alpha_A_PC1_R <- extract_b_parameter(d_pc1, "X_adjacent")
  alpha_F_PC1_R <- extract_b_parameter(d_pc1, "X_far")
  alpha_A_PC2_R <- extract_b_parameter(d_pc2, "X_adjacent")
  alpha_F_PC2_R <- extract_b_parameter(d_pc2, "X_far")
  alpha_A_PC3 <- extract_b_parameter(d_pc3, "X_adjacent")
  alpha_F_PC3 <- extract_b_parameter(d_pc3, "X_far")

  beta_PC1_R <- extract_b_parameter(d_out, "PC1_R")
  beta_PC2_R <- extract_b_parameter(d_out, "PC2_R")
  beta_PC3 <- extract_b_parameter(d_out, "PC3")

  IE_PC1_A <- alpha_A_PC1_R * beta_PC1_R
  IE_PC2_A <- alpha_A_PC2_R * beta_PC2_R
  IE_PC3_A <- alpha_A_PC3 * beta_PC3
  IE_PC1_F <- alpha_F_PC1_R * beta_PC1_R
  IE_PC2_F <- alpha_F_PC2_R * beta_PC2_R
  IE_PC3_F <- alpha_F_PC3 * beta_PC3

  NIE_A <- IE_PC1_A + IE_PC2_A + IE_PC3_A
  NIE_F <- IE_PC1_F + IE_PC2_F + IE_PC3_F

  DE_A <- extract_b_parameter(d_out, "X_adjacent")
  DE_F <- extract_b_parameter(d_out, "X_far")
  TE_A <- extract_b_parameter(d_total, "X_adjacent")
  TE_F <- extract_b_parameter(d_total, "X_far")

  quantity_draws <- tibble::tibble(
    draw = seq_len(n_draws),
    alpha_A_PC1_R = alpha_A_PC1_R,
    alpha_F_PC1_R = alpha_F_PC1_R,
    alpha_A_PC2_R = alpha_A_PC2_R,
    alpha_F_PC2_R = alpha_F_PC2_R,
    alpha_A_PC3 = alpha_A_PC3,
    alpha_F_PC3 = alpha_F_PC3,
    beta_PC1_R = beta_PC1_R,
    beta_PC2_R = beta_PC2_R,
    beta_PC3 = beta_PC3,
    IE_PC1_A = IE_PC1_A,
    IE_PC2_A = IE_PC2_A,
    IE_PC3_A = IE_PC3_A,
    IE_PC1_F = IE_PC1_F,
    IE_PC2_F = IE_PC2_F,
    IE_PC3_F = IE_PC3_F,
    NIE_A = NIE_A,
    NIE_F = NIE_F,
    DE_A = DE_A,
    DE_F = DE_F,
    TE_A = TE_A,
    TE_F = TE_F,
    PM_A = NIE_A / TE_A,
    PM_F = NIE_F / TE_F
  )

  list(
    draws = quantity_draws,
    sampled_indices = idx,
    available_draws = available_draws,
    n_draws = n_draws,
    path_summary = summarize_bayesian_paths(list(draws = quantity_draws)),
    decomposition_summary = summarize_bayesian_decomposition(list(draws = quantity_draws))
  )
}

compute_frequentist_current_estimates <- function(data,
                                                  exposures = c("X_adjacent", "X_far"),
                                                  outcome = "Y",
                                                  mediators = c("PC1_R", "PC2_R", "PC3"),
                                                  covariates = c("x_coord", "y_coord")) {
  if (!exists("decompose_linear_multix", mode = "function")) {
    stop("decompose_linear_multix() is missing. Source R/pca_mediation_pipeline.R first.", call. = FALSE)
  }

  fit <- decompose_linear_multix(
    data = data,
    exposures = exposures,
    outcome = outcome,
    mediators = mediators,
    covariates = covariates
  )

  list(
    fit = fit,
    summary = fit$summary |>
      dplyr::mutate(
        contrast = dplyr::case_when(
          exposure == "X_adjacent" ~ "adjacent vs inside",
          exposure == "X_far" ~ "far vs inside",
          TRUE ~ exposure
        )
      ),
    path = fit$path |>
      dplyr::mutate(
        contrast = dplyr::case_when(
          exposure == "X_adjacent" ~ "adjacent vs inside",
          exposure == "X_far" ~ "far vs inside",
          TRUE ~ exposure
        )
      )
  )
}

compare_frequentist_bayesian_v0 <- function(frequentist,
                                            bayesian_path_summary,
                                            bayesian_decomposition_summary) {
  decomposition_comparison <- frequentist$summary |>
    dplyr::select(contrast, frequentist_TE = TE, frequentist_direct = NDE, frequentist_NIE = NIE) |>
    dplyr::left_join(
      bayesian_decomposition_summary |>
        dplyr::select(
          contrast,
          bayesian_TE_mean = TE_mean,
          bayesian_direct_mean = direct_mean,
          bayesian_NIE_mean = NIE_mean
        ),
      by = "contrast"
    )

  path_comparison <- frequentist$path |>
    dplyr::select(
      contrast,
      mediator,
      frequentist_alpha = alpha_X_to_M,
      frequentist_beta = beta_M_to_Y,
      frequentist_indirect = indirect_component
    ) |>
    dplyr::left_join(
      bayesian_path_summary |>
        dplyr::select(
          contrast,
          mediator,
          bayesian_alpha_mean = alpha_mean,
          bayesian_beta_mean = beta_mean,
          bayesian_indirect_mean = indirect_mean
        ),
      by = c("contrast", "mediator")
    )

  list(
    decomposition = decomposition_comparison,
    path = path_comparison
  )
}

summarize_bayesian_diagnostics <- function(fits, max_treedepth = 10) {
  check_bayesian_v0_packages()

  purrr::imap_dfr(fits, function(fit, model) {
    draw_diag <- posterior::summarise_draws(
      posterior::as_draws(fit),
      "rhat",
      "ess_bulk",
      "ess_tail"
    )

    nuts <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
    divergent_transitions <- NA_integer_
    max_treedepth_hits <- NA_integer_

    if (!is.null(nuts)) {
      divergent_transitions <- sum(nuts$Parameter == "divergent__" & nuts$Value == 1, na.rm = TRUE)
      max_treedepth_hits <- sum(nuts$Parameter == "treedepth__" & nuts$Value >= max_treedepth, na.rm = TRUE)
    }

    tibble::tibble(
      model = model,
      max_Rhat = max(draw_diag$rhat, na.rm = TRUE),
      min_bulk_ESS = min(draw_diag$ess_bulk, na.rm = TRUE),
      min_tail_ESS = min(draw_diag$ess_tail, na.rm = TRUE),
      divergent_transitions = divergent_transitions,
      max_treedepth_hits = max_treedepth_hits,
      max_treedepth_issue = !is.na(max_treedepth_hits) && max_treedepth_hits > 0
    )
  })
}
