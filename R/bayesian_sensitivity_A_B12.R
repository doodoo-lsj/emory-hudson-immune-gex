# Bayesian sensitivity analyses A, B1, and B2 for the current PCA mediation model.
#
# This file intentionally excludes latent-U B3 and B4. Mediator-outcome
# confounding is handled only by Sensitivity A, using fixed residual
# correlations rho_j.

check_bayesian_sensitivity_packages <- function(require_brms = TRUE,
                                                require_stan = FALSE) {
  required <- c("posterior", "dplyr", "tidyr", "tibble", "purrr")
  if (isTRUE(require_brms)) {
    required <- c("brms", required)
  }
  if (isTRUE(require_stan)) {
    required <- c("rstan", required)
  }
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Bayesian sensitivity requires installed R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

bayes_sens_contrast_label <- function(exposure) {
  dplyr::case_when(
    exposure == "X_adjacent" ~ "adjacent vs inside",
    exposure == "X_far" ~ "far vs inside",
    TRUE ~ exposure
  )
}

bayes_sens_coord_terms <- function(standardize_covariates = TRUE) {
  if (isTRUE(standardize_covariates)) c("x_std", "y_std") else c("x_coord", "y_coord")
}

bayes_sens_summary <- function(x) {
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

bayes_sens_summarize_draws <- function(draws, value_col = "value", group_cols) {
  draws |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      mean = mean(.data[[value_col]], na.rm = TRUE),
      median = stats::median(.data[[value_col]], na.rm = TRUE),
      sd = stats::sd(.data[[value_col]], na.rm = TRUE),
      q025 = stats::quantile(.data[[value_col]], 0.025, na.rm = TRUE),
      q975 = stats::quantile(.data[[value_col]], 0.975, na.rm = TRUE),
      Pr_gt_0 = mean(.data[[value_col]] > 0, na.rm = TRUE),
      Pr_lt_0 = mean(.data[[value_col]] < 0, na.rm = TRUE),
      .groups = "drop"
    )
}

bayes_sens_extract_b <- function(draws_df, coef_name, resp = NULL) {
  extract_bayes_v1_b(draws_df, coef_name = coef_name, resp = resp)
}

bayes_sens_extract_sigma <- function(draws_df, resp) {
  extract_bayes_v1_sigma(draws_df, resp = resp)
}

bayes_sens_extract_univariate_sigma <- function(draws_df, label = "univariate outcome") {
  candidates <- c("sigma", "sigma_Y", "sigma.Y")
  existing <- intersect(candidates, names(draws_df))
  if (length(existing) == 0) {
    stop("Could not find posterior sigma column for ", label, ". Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  draws_df[[existing[[1]]]]
}

bayes_sens_extract_rescor <- function(draws_df, resp1, resp2) {
  extract_bayes_v1_rescor(draws_df, resp1 = resp1, resp2 = resp2)
}

compute_bayesian_baseline_decomposition_draws <- function(fits = NULL,
                                                          n_draws = NULL,
                                                          seed = 123,
                                                          residual_draws = NULL) {
  check_bayesian_sensitivity_packages(require_brms = TRUE)
  if (is.null(residual_draws)) {
    residual_draws <- extract_bayesian_A_residual_draws(
      fits = fits,
      n_draws = n_draws,
      seed = seed
    )
  }
  n_draws <- residual_draws$n_draws
  alpha <- residual_draws$alpha
  beta <- residual_draws$beta
  direct <- residual_draws$direct

  make_wide <- function(exposure) {
    alpha_mat <- alpha[[exposure]]
    direct_vec <- direct[, exposure]
    IIE <- alpha_mat * beta
    total_IIE <- rowSums(IIE)
    TE <- direct_vec + total_IIE
    PM <- total_IIE / TE
    tibble::tibble(
      draw = seq_len(n_draws),
      exposure = exposure,
      contrast = bayes_sens_contrast_label(exposure),
      direct_effect = direct_vec,
      total_IIE = total_IIE,
      TE = TE,
      PM = PM,
      IIE_PC1_R = IIE[, "PC1_R"],
      IIE_PC2_R = IIE[, "PC2_R"],
      IIE_PC3 = IIE[, "PC3"],
      alpha_PC1_R = alpha_mat[, "PC1_R"],
      alpha_PC2_R = alpha_mat[, "PC2_R"],
      alpha_PC3 = alpha_mat[, "PC3"],
      beta_PC1_R = beta[, "PC1_R"],
      beta_PC2_R = beta[, "PC2_R"],
      beta_PC3 = beta[, "PC3"]
    )
  }
  wide <- dplyr::bind_rows(make_wide("X_adjacent"), make_wide("X_far"))

  path_draws <- wide |>
    dplyr::select(draw, exposure, contrast, dplyr::starts_with("alpha_"), dplyr::starts_with("beta_"), dplyr::starts_with("IIE_")) |>
    tidyr::pivot_longer(
      cols = dplyr::matches("^(alpha|beta|IIE)_"),
      names_to = c(".value", "mediator"),
      names_pattern = "^(alpha|beta|IIE)_(.*)$"
    )

  decomposition_draws <- wide |>
    dplyr::select(draw, exposure, contrast, direct_effect, total_IIE, TE, PM)

  list(
    wide = wide,
    path_draws = path_draws,
    decomposition_draws = decomposition_draws,
    path_summary = path_draws |>
      dplyr::group_by(contrast, exposure, mediator) |>
      dplyr::summarise(
        alpha_mean = mean(alpha), alpha_median = stats::median(alpha),
        alpha_q025 = stats::quantile(alpha, 0.025), alpha_q975 = stats::quantile(alpha, 0.975),
        beta_mean = mean(beta), beta_median = stats::median(beta),
        beta_q025 = stats::quantile(beta, 0.025), beta_q975 = stats::quantile(beta, 0.975),
        IIE_mean = mean(IIE), IIE_median = stats::median(IIE),
        IIE_q025 = stats::quantile(IIE, 0.025), IIE_q975 = stats::quantile(IIE, 0.975),
        Pr_IIE_gt_0 = mean(IIE > 0),
        .groups = "drop"
      ),
    decomposition_summary = decomposition_draws |>
      tidyr::pivot_longer(cols = c(direct_effect, total_IIE, TE, PM), names_to = "quantity", values_to = "value") |>
      bayes_sens_summarize_draws(group_cols = c("contrast", "exposure", "quantity")) |>
      tidyr::pivot_wider(names_from = quantity, values_from = c(mean, median, sd, q025, q975, Pr_gt_0, Pr_lt_0), names_glue = "{quantity}_{.value}"),
    residual_draws = residual_draws,
    sampled_indices = residual_draws$sampled_indices,
    n_draws = n_draws
  )
}

extract_bayesian_A_residual_draws <- function(fits,
                                              n_draws = NULL,
                                              seed = 123,
                                              sampled_indices = NULL) {
  d_m_all <- posterior::as_draws_df(fits$mediator_joint)
  d_y_all <- posterior::as_draws_df(fits$outcome)
  available <- c(mediator_joint = nrow(d_m_all), outcome = nrow(d_y_all))
  if (is.null(sampled_indices)) {
    set.seed(seed)
    if (is.null(n_draws)) n_draws <- min(available)
    idx_m <- sample.int(nrow(d_m_all), n_draws, replace = TRUE)
    idx_y <- sample.int(nrow(d_y_all), n_draws, replace = TRUE)
  } else {
    idx_m <- sampled_indices$mediator_joint
    idx_y <- sampled_indices$outcome
    n_draws <- length(idx_m)
    if (length(idx_y) != n_draws) {
      stop("Mediator and outcome sampled index vectors must have the same length.", call. = FALSE)
    }
  }
  d_m <- d_m_all[idx_m, , drop = FALSE]
  d_y <- d_y_all[idx_y, , drop = FALSE]

  sig <- cbind(
    PC1_R = bayes_sens_extract_sigma(d_m, "PC1R"),
    PC2_R = bayes_sens_extract_sigma(d_m, "PC2R"),
    PC3 = bayes_sens_extract_sigma(d_m, "PC3")
  )
  r12 <- bayes_sens_extract_rescor(d_m, "PC1R", "PC2R")
  r13 <- bayes_sens_extract_rescor(d_m, "PC1R", "PC3")
  r23 <- bayes_sens_extract_rescor(d_m, "PC2R", "PC3")
  sigma_y <- bayes_sens_extract_univariate_sigma(d_y, "outcome model")

  beta <- cbind(
    PC1_R = bayes_sens_extract_b(d_y, "PC1_R"),
    PC2_R = bayes_sens_extract_b(d_y, "PC2_R"),
    PC3 = bayes_sens_extract_b(d_y, "PC3")
  )
  alpha_A <- cbind(
    PC1_R = bayes_sens_extract_b(d_m, "X_adjacent", "PC1R"),
    PC2_R = bayes_sens_extract_b(d_m, "X_adjacent", "PC2R"),
    PC3 = bayes_sens_extract_b(d_m, "X_adjacent", "PC3")
  )
  alpha_F <- cbind(
    PC1_R = bayes_sens_extract_b(d_m, "X_far", "PC1R"),
    PC2_R = bayes_sens_extract_b(d_m, "X_far", "PC2R"),
    PC3 = bayes_sens_extract_b(d_m, "X_far", "PC3")
  )
  direct <- cbind(
    X_adjacent = bayes_sens_extract_b(d_y, "X_adjacent"),
    X_far = bayes_sens_extract_b(d_y, "X_far")
  )

  list(
    n_draws = n_draws,
    sig = sig,
    r12 = r12,
    r13 = r13,
    r23 = r23,
    sigma_y = sigma_y,
    beta = beta,
    alpha = list(X_adjacent = alpha_A, X_far = alpha_F),
    direct = direct,
    sampled_indices = list(mediator_joint = idx_m, outcome = idx_y)
  )
}

bayesian_A_rho_grid <- function(rho_grid = seq(-0.5, 0.5, by = 0.05),
                                mode = c("common", "one_at_a_time"),
                                mediators = c("PC1_R", "PC2_R", "PC3")) {
  mode <- match.arg(mode)
  if (mode == "common") {
    return(tibble::tibble(scenario = "common", mediator_varied = "all", rho = rho_grid, rho_1 = rho_grid, rho_2 = rho_grid, rho_3 = rho_grid))
  }
  purrr::map_dfr(mediators, function(m) {
    idx <- match(m, mediators)
    out <- tibble::tibble(scenario = "one_at_a_time", mediator_varied = m, rho = rho_grid, rho_1 = 0, rho_2 = 0, rho_3 = 0)
    out[[paste0("rho_", idx)]] <- rho_grid
    out
  })
}

compute_bayesian_sensitivity_A <- function(fits,
                                           baseline = NULL,
                                           residual_draws = NULL,
                                           rho_grid = seq(-0.5, 0.5, by = 0.05),
                                           mode = c("common", "one_at_a_time"),
                                           n_draws = NULL,
                                           seed = 123,
                                           feasibility_tol = 1e-10) {
  mode <- match.arg(mode)
  if (is.null(residual_draws) && !is.null(baseline$residual_draws)) {
    residual_draws <- baseline$residual_draws
  }
  if (is.null(residual_draws)) {
    residual_draws <- extract_bayesian_A_residual_draws(fits, n_draws = n_draws, seed = seed)
  }
  grid <- bayesian_A_rho_grid(rho_grid, mode = mode)
  mediator_names <- c("PC1_R", "PC2_R", "PC3")

  rows <- purrr::pmap_dfr(grid, function(scenario, mediator_varied, rho, rho_1, rho_2, rho_3) {
    rho_vec <- c(rho_1, rho_2, rho_3)
    draw_rows <- vector("list", residual_draws$n_draws)
    for (s in seq_len(residual_draws$n_draws)) {
      Sigma_M <- diag(residual_draws$sig[s, ]^2, 3)
      Sigma_M[1, 2] <- Sigma_M[2, 1] <- residual_draws$r12[[s]] * residual_draws$sig[s, 1] * residual_draws$sig[s, 2]
      Sigma_M[1, 3] <- Sigma_M[3, 1] <- residual_draws$r13[[s]] * residual_draws$sig[s, 1] * residual_draws$sig[s, 3]
      Sigma_M[2, 3] <- Sigma_M[3, 2] <- residual_draws$r23[[s]] * residual_draws$sig[s, 2] * residual_draws$sig[s, 3]
      c_my <- rho_vec * residual_draws$sig[s, ] * residual_draws$sigma_y[[s]]
      feasibility <- as.numeric(t(c_my) %*% solve(Sigma_M, c_my) / residual_draws$sigma_y[[s]]^2)
      valid <- is.finite(feasibility) && feasibility < (1 - feasibility_tol)
      beta_adj <- rep(NA_real_, 3)
      if (valid) {
        beta_adj <- as.numeric(residual_draws$beta[s, ] - solve(Sigma_M, c_my))
      }
      names(beta_adj) <- mediator_names
      draw_rows[[s]] <- tibble::tibble(
        draw = s,
        valid_parameter = valid,
        feasibility = feasibility,
        beta_PC1_R = beta_adj[[1]],
        beta_PC2_R = beta_adj[[2]],
        beta_PC3 = beta_adj[[3]]
      )
    }
    dplyr::bind_rows(draw_rows) |>
      dplyr::mutate(scenario = scenario, mediator_varied = mediator_varied, rho = rho, rho_1 = rho_1, rho_2 = rho_2, rho_3 = rho_3, .before = draw)
  })

  make_contrast <- function(exposure) {
    alpha <- residual_draws$alpha[[exposure]]
    direct <- residual_draws$direct[, exposure]
    rows |>
      dplyr::mutate(
        exposure = exposure,
        contrast = bayes_sens_contrast_label(exposure),
        IIE_PC1_R = alpha[draw, "PC1_R"] * beta_PC1_R,
        IIE_PC2_R = alpha[draw, "PC2_R"] * beta_PC2_R,
        IIE_PC3 = alpha[draw, "PC3"] * beta_PC3,
        total_IIE = IIE_PC1_R + IIE_PC2_R + IIE_PC3,
        direct_effect = direct[draw],
        TE = direct_effect + total_IIE,
        PM = total_IIE / TE,
        alpha_PC1_R = alpha[draw, "PC1_R"],
        alpha_PC2_R = alpha[draw, "PC2_R"],
        alpha_PC3 = alpha[draw, "PC3"]
      )
  }
  wide <- dplyr::bind_rows(make_contrast("X_adjacent"), make_contrast("X_far"))

  path_draws <- wide |>
    dplyr::select(scenario, mediator_varied, rho, rho_1, rho_2, rho_3, draw, valid_parameter, exposure, contrast, dplyr::starts_with("alpha_"), dplyr::starts_with("beta_"), dplyr::starts_with("IIE_")) |>
    tidyr::pivot_longer(
      cols = dplyr::matches("^(alpha|beta|IIE)_"),
      names_to = c(".value", "mediator"),
      names_pattern = "^(alpha|beta|IIE)_(.*)$"
    )
  summary_draws <- wide |>
    dplyr::select(scenario, mediator_varied, rho, rho_1, rho_2, rho_3, draw, valid_parameter, exposure, contrast, direct_effect, total_IIE, TE, PM)

  list(
    grid = grid,
    residual_draws = residual_draws,
    draws_wide = wide,
    path = path_draws |>
      dplyr::filter(valid_parameter) |>
      dplyr::group_by(scenario, mediator_varied, rho, rho_1, rho_2, rho_3, contrast, exposure, mediator) |>
      dplyr::summarise(
        alpha_mean = mean(alpha), alpha_q025 = stats::quantile(alpha, 0.025), alpha_q975 = stats::quantile(alpha, 0.975),
        beta_mean = mean(beta), beta_q025 = stats::quantile(beta, 0.025), beta_q975 = stats::quantile(beta, 0.975),
        IIE_mean = mean(IIE), IIE_median = stats::median(IIE), IIE_q025 = stats::quantile(IIE, 0.025), IIE_q975 = stats::quantile(IIE, 0.975),
        .groups = "drop"
      ),
    summary = summary_draws |>
      dplyr::filter(valid_parameter) |>
      tidyr::pivot_longer(cols = c(direct_effect, total_IIE, TE, PM), names_to = "quantity", values_to = "value") |>
      bayes_sens_summarize_draws(group_cols = c("scenario", "mediator_varied", "rho", "rho_1", "rho_2", "rho_3", "contrast", "exposure", "quantity")) |>
      tidyr::pivot_wider(names_from = quantity, values_from = c(mean, median, sd, q025, q975, Pr_gt_0, Pr_lt_0), names_glue = "{quantity}_{.value}"),
    feasibility = wide |>
      dplyr::group_by(scenario, mediator_varied, rho, rho_1, rho_2, rho_3) |>
      dplyr::summarise(valid_draw_fraction = mean(valid_parameter), max_feasibility = max(feasibility, na.rm = TRUE), .groups = "drop"),
    validation = list(
      TE_identity_max_abs = max(abs(wide$TE - wide$direct_effect - wide$total_IIE), na.rm = TRUE),
      rho0_available = any(grid$rho_1 == 0 & grid$rho_2 == 0 & grid$rho_3 == 0)
    )
  )
}

compute_bayesian_sensitivity_A_suite <- function(fits,
                                                 n_draws = NULL,
                                                 seed = 123,
                                                 rho_grid = c(-0.3, 0, 0.3),
                                                 feasibility_tol = 1e-10) {
  shared_A_draws <- extract_bayesian_A_residual_draws(
    fits = fits,
    n_draws = n_draws,
    seed = seed
  )
  baseline_quantities <- compute_bayesian_baseline_decomposition_draws(
    residual_draws = shared_A_draws
  )
  A_common <- compute_bayesian_sensitivity_A(
    fits = fits,
    baseline = baseline_quantities,
    residual_draws = shared_A_draws,
    rho_grid = rho_grid,
    mode = "common",
    feasibility_tol = feasibility_tol
  )
  A_one <- compute_bayesian_sensitivity_A(
    fits = fits,
    baseline = baseline_quantities,
    residual_draws = shared_A_draws,
    rho_grid = rho_grid,
    mode = "one_at_a_time",
    feasibility_tol = feasibility_tol
  )
  validation <- list(
    common_TE_identity = A_common$validation$TE_identity_max_abs,
    one_at_a_time_TE_identity = A_one$validation$TE_identity_max_abs,
    common_rho0_vs_baseline = validate_bayesian_A_rho0_against_baseline(A_common, baseline_quantities),
    one_at_a_time_rho0_vs_baseline = validate_bayesian_A_rho0_against_baseline(A_one, baseline_quantities),
    all_A_rho0_scenarios_identical = validate_bayesian_A_rho0_scenarios_identical(A_common, A_one),
    sampled_indices_identical = list(
      baseline_common = identical(baseline_quantities$sampled_indices, A_common$residual_draws$sampled_indices),
      baseline_one_at_a_time = identical(baseline_quantities$sampled_indices, A_one$residual_draws$sampled_indices),
      common_one_at_a_time = identical(A_common$residual_draws$sampled_indices, A_one$residual_draws$sampled_indices)
    ),
    common_feasibility = A_common$feasibility,
    one_at_a_time_feasibility = A_one$feasibility
  )
  list(
    baseline_quantities = baseline_quantities,
    common = A_common,
    one_at_a_time = A_one,
    validation = validation,
    shared_draws = shared_A_draws
  )
}

bayes_sens_B_stan_code <- function() {
  "
data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> Pm;
  int<lower=1> Py;
  matrix[N, Pm] X_m;
  matrix[N, Py] X_y;
  matrix[N, K] M;
  vector[N] Y;
  vector<lower=0, upper=1>[N] pi_u;
  vector[K] lambda_U;
  real delta_U;
  matrix[Pm, K] prior_loc_B_m;
  matrix<lower=0>[Pm, K] prior_scale_B_m;
  vector[Py] prior_loc_b_y;
  vector<lower=0>[Py] prior_scale_b_y;
  vector<lower=0>[K] rate_sigma_m;
  real<lower=0> rate_sigma_y;
}
parameters {
  matrix[Pm, K] B_m;
  vector[Py] b_y;
  vector<lower=0>[K] sigma_m;
  real<lower=0> sigma_y;
  cholesky_factor_corr[K] Lcorr_m;
}
transformed parameters {
  matrix[K, K] L_Sigma_m;
  L_Sigma_m = diag_pre_multiply(sigma_m, Lcorr_m);
}
model {
  for (k in 1:K) {
    for (p in 1:Pm) {
      B_m[p, k] ~ student_t(3, prior_loc_B_m[p, k], prior_scale_B_m[p, k]);
    }
    sigma_m[k] ~ exponential(rate_sigma_m[k]);
  }
  for (p in 1:Py) {
    b_y[p] ~ student_t(3, prior_loc_b_y[p], prior_scale_b_y[p]);
  }
  sigma_y ~ exponential(rate_sigma_y);
  Lcorr_m ~ lkj_corr_cholesky(2);

  for (n in 1:N) {
    row_vector[K] mu_m0_row = X_m[n] * B_m;
    vector[K] mu_m0 = to_vector(mu_m0_row);
    vector[K] mu_m1 = mu_m0 - lambda_U;
    real mu_y0 = dot_product(to_vector(X_y[n]), b_y);
    real mu_y1 = mu_y0 - delta_U;
    real ll0 = log1m(pi_u[n]) +
      multi_normal_cholesky_lpdf(to_vector(M[n]) | mu_m0, L_Sigma_m) +
      normal_lpdf(Y[n] | mu_y0, sigma_y);
    real ll1 = log(pi_u[n]) +
      multi_normal_cholesky_lpdf(to_vector(M[n]) | mu_m1, L_Sigma_m) +
      normal_lpdf(Y[n] | mu_y1, sigma_y);
    target += log_sum_exp(ll0, ll1);
  }
}
"
}

.bayes_sens_stan_cache <- new.env(parent = emptyenv())

bayes_sens_compile_B_model <- function(force_recompile = FALSE) {
  check_bayesian_sensitivity_packages(require_brms = FALSE, require_stan = TRUE)
  if (!force_recompile && exists("B_model", envir = .bayes_sens_stan_cache, inherits = FALSE)) {
    return(get("B_model", envir = .bayes_sens_stan_cache))
  }
  rstan::rstan_options(auto_write = TRUE)
  model <- rstan::stan_model(model_code = bayes_sens_B_stan_code(), model_name = "bayes_sensitivity_B12")
  assign("B_model", model, envir = .bayes_sens_stan_cache)
  model
}

bayes_sens_prior_data <- function(data, formula, response) {
  y <- data[[response]]
  y_sd <- stats::sd(y)
  y_mean <- mean(y)
  x <- stats::model.matrix(formula, data = data)
  loc <- rep(0, ncol(x))
  loc[colnames(x) == "(Intercept)"] <- y_mean
  scale <- vapply(seq_len(ncol(x)), function(j) {
    if (colnames(x)[[j]] == "(Intercept)") {
      10 * y_sd
    } else {
      x_sd <- stats::sd(x[, j])
      10 * y_sd / x_sd
    }
  }, numeric(1))
  list(loc = loc, scale = scale, y_sd = y_sd, y_mean = y_mean, x = x)
}

prepare_bayesian_B_stan_data <- function(data,
                                         scenario = c("B1_XY", "B2_XM"),
                                         p_U = 0.3,
                                         OR_U_adj = 1,
                                         OR_U_far = 1,
                                         s_M = c(0, 0, 0),
                                         s_Y = 0,
                                         standardize_covariates = TRUE) {
  scenario <- match.arg(scenario)
  coord_terms <- bayes_sens_coord_terms(standardize_covariates)
  m_formula <- stats::as.formula(paste("PC1_R ~ X_adjacent + X_far +", paste(coord_terms, collapse = " + ")))
  y_formula <- stats::as.formula(paste("Y ~ X_adjacent + X_far + PC1_R + PC2_R + PC3 +", paste(coord_terms, collapse = " + ")))
  X_m <- stats::model.matrix(m_formula, data = data)
  X_y <- stats::model.matrix(y_formula, data = data)
  M <- as.matrix(data[, c("PC1_R", "PC2_R", "PC3")])
  Y <- data$Y

  fixed <- sensitivity_B_make_fixed_params(
    data = data,
    p_U = p_U,
    OR_U_adj = OR_U_adj,
    OR_U_far = OR_U_far,
    s_M = if (scenario == "B1_XY") c(0, 0, 0) else s_M,
    s_Y = if (scenario == "B2_XM") 0 else s_Y
  )
  pi_u <- sensitivity_B_pi(data, fixed)

  prior_m <- lapply(c("PC1_R", "PC2_R", "PC3"), function(resp) bayes_sens_prior_data(data, m_formula, resp))
  prior_y <- bayes_sens_prior_data(data, y_formula, "Y")
  prior_loc_B_m <- do.call(cbind, lapply(prior_m, `[[`, "loc"))
  prior_scale_B_m <- do.call(cbind, lapply(prior_m, `[[`, "scale"))
  prior_loc_b_y <- prior_y$loc
  prior_scale_b_y <- prior_y$scale

  list(
    stan_data = list(
      N = nrow(data),
      K = 3,
      Pm = ncol(X_m),
      Py = ncol(X_y),
      X_m = X_m,
      X_y = X_y,
      M = M,
      Y = Y,
      pi_u = as.numeric(pi_u),
      lambda_U = as.numeric(fixed$lambda_U),
      delta_U = fixed$delta_U,
      prior_loc_B_m = prior_loc_B_m,
      prior_scale_B_m = prior_scale_B_m,
      prior_loc_b_y = prior_loc_b_y,
      prior_scale_b_y = prior_scale_b_y,
      rate_sigma_m = 1 / vapply(prior_m, `[[`, numeric(1), "y_sd"),
      rate_sigma_y = 1 / prior_y$y_sd
    ),
    fixed = fixed,
    terms = list(
      mediator = colnames(X_m),
      outcome = colnames(X_y),
      mediators = c("PC1_R", "PC2_R", "PC3")
    ),
    scenario = scenario
  )
}

fit_bayesian_B_scenario <- function(data,
                                    scenario = c("B1_XY", "B2_XM"),
                                    p_U = 0.3,
                                    OR_U_adj = 1,
                                    OR_U_far = 1,
                                    s_M = c(0, 0, 0),
                                    s_Y = 0,
                                    standardize_covariates = TRUE,
                                    chains = 4,
                                    iter = 1000,
                                    warmup = 500,
                                    cores = 4,
                                    seed = 123,
                                    adapt_delta = 0.9,
                                    max_treedepth = 10,
                                    refresh = 100,
                                    stan_model = NULL) {
  scenario <- match.arg(scenario)
  check_bayesian_sensitivity_packages(require_brms = FALSE, require_stan = TRUE)
  prepared <- prepare_bayesian_B_stan_data(
    data = data,
    scenario = scenario,
    p_U = p_U,
    OR_U_adj = OR_U_adj,
    OR_U_far = OR_U_far,
    s_M = s_M,
    s_Y = s_Y,
    standardize_covariates = standardize_covariates
  )
  if (is.null(stan_model)) stan_model <- bayes_sens_compile_B_model()
  elapsed <- system.time({
    fit <- rstan::sampling(
      object = stan_model,
      data = prepared$stan_data,
      chains = chains,
      iter = iter,
      warmup = warmup,
      cores = cores,
      seed = seed,
      refresh = refresh,
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
    )
  })
  list(fit = fit, prepared = prepared, timing = tibble::tibble(component = paste0("fit_", scenario), elapsed_seconds = unname(elapsed[["elapsed"]])))
}

extract_rstan_draw <- function(draws_df, name) {
  if (!name %in% names(draws_df)) {
    stop("Missing Stan draw column: ", name, call. = FALSE)
  }
  draws_df[[name]]
}

compute_bayesian_B_quantities <- function(fit_object, n_draws = NULL, seed = 123) {
  set.seed(seed)
  fit <- fit_object$fit
  prepared <- fit_object$prepared
  draws_all <- posterior::as_draws_df(fit)
  if (is.null(n_draws)) n_draws <- nrow(draws_all)
  idx <- sample.int(nrow(draws_all), n_draws, replace = TRUE)
  d <- draws_all[idx, , drop = FALSE]
  terms_m <- prepared$terms$mediator
  terms_y <- prepared$terms$outcome
  meds <- prepared$terms$mediators

  get_B <- function(term, med) {
    p <- match(term, terms_m)
    k <- match(med, meds)
    extract_rstan_draw(d, paste0("B_m[", p, ",", k, "]"))
  }
  get_y <- function(term) {
    p <- match(term, terms_y)
    extract_rstan_draw(d, paste0("b_y[", p, "]"))
  }

  beta <- cbind(PC1_R = get_y("PC1_R"), PC2_R = get_y("PC2_R"), PC3 = get_y("PC3"))
  direct <- cbind(X_adjacent = get_y("X_adjacent"), X_far = get_y("X_far"))

  make_contrast <- function(exposure) {
    alpha <- cbind(
      PC1_R = get_B(exposure, "PC1_R"),
      PC2_R = get_B(exposure, "PC2_R"),
      PC3 = get_B(exposure, "PC3")
    )
    IIE <- alpha * beta
    total_IIE <- rowSums(IIE)
    TE <- direct[, exposure] + total_IIE
    tibble::tibble(
      draw = seq_len(n_draws),
      scenario = prepared$scenario,
      p_U = prepared$fixed$p_U,
      OR_U_adj = prepared$fixed$OR_U_adj,
      OR_U_far = prepared$fixed$OR_U_far,
      s_M1 = prepared$fixed$s_M[[1]],
      s_M2 = prepared$fixed$s_M[[2]],
      s_M3 = prepared$fixed$s_M[[3]],
      s_Y = prepared$fixed$s_Y,
      exposure = exposure,
      contrast = bayes_sens_contrast_label(exposure),
      direct_effect = direct[, exposure],
      total_IIE = total_IIE,
      TE = TE,
      PM = total_IIE / TE,
      alpha_PC1_R = alpha[, "PC1_R"], alpha_PC2_R = alpha[, "PC2_R"], alpha_PC3 = alpha[, "PC3"],
      beta_PC1_R = beta[, "PC1_R"], beta_PC2_R = beta[, "PC2_R"], beta_PC3 = beta[, "PC3"],
      IIE_PC1_R = IIE[, "PC1_R"], IIE_PC2_R = IIE[, "PC2_R"], IIE_PC3 = IIE[, "PC3"]
    )
  }
  wide <- dplyr::bind_rows(make_contrast("X_adjacent"), make_contrast("X_far"))

  path_draws <- wide |>
    dplyr::select(scenario, p_U, OR_U_adj, OR_U_far, s_M1, s_M2, s_M3, s_Y, draw, exposure, contrast, dplyr::starts_with("alpha_"), dplyr::starts_with("beta_"), dplyr::starts_with("IIE_")) |>
    tidyr::pivot_longer(cols = dplyr::matches("^(alpha|beta|IIE)_"), names_to = c(".value", "mediator"), names_pattern = "^(alpha|beta|IIE)_(.*)$")
  summary_draws <- wide |>
    dplyr::select(scenario, p_U, OR_U_adj, OR_U_far, s_M1, s_M2, s_M3, s_Y, draw, exposure, contrast, direct_effect, total_IIE, TE, PM)

  list(
    draws_wide = wide,
    path_draws = path_draws,
    decomposition_draws = summary_draws,
    path = path_draws |>
      dplyr::group_by(scenario, p_U, OR_U_adj, OR_U_far, s_M1, s_M2, s_M3, s_Y, contrast, exposure, mediator) |>
      dplyr::summarise(
        alpha_mean = mean(alpha), alpha_q025 = stats::quantile(alpha, 0.025), alpha_q975 = stats::quantile(alpha, 0.975),
        beta_mean = mean(beta), beta_q025 = stats::quantile(beta, 0.025), beta_q975 = stats::quantile(beta, 0.975),
        IIE_mean = mean(IIE), IIE_median = stats::median(IIE), IIE_q025 = stats::quantile(IIE, 0.025), IIE_q975 = stats::quantile(IIE, 0.975),
        .groups = "drop"
      ),
    summary = summary_draws |>
      tidyr::pivot_longer(cols = c(direct_effect, total_IIE, TE, PM), names_to = "quantity", values_to = "value") |>
      bayes_sens_summarize_draws(group_cols = c("scenario", "p_U", "OR_U_adj", "OR_U_far", "s_M1", "s_M2", "s_M3", "s_Y", "contrast", "exposure", "quantity")) |>
      tidyr::pivot_wider(names_from = quantity, values_from = c(mean, median, sd, q025, q975, Pr_gt_0, Pr_lt_0), names_glue = "{quantity}_{.value}"),
    validation = list(TE_identity_max_abs = max(abs(wide$TE - wide$direct_effect - wide$total_IIE), na.rm = TRUE)),
    sampled_indices = idx,
    n_draws = n_draws
  )
}

summarize_rstan_sensitivity_diagnostics <- function(fit, max_treedepth = 10) {
  draw_diag <- posterior::summarise_draws(posterior::as_draws(fit), "rhat", "ess_bulk", "ess_tail")
  sampler <- tryCatch(rstan::get_sampler_params(fit, inc_warmup = FALSE), error = function(e) NULL)
  divergences <- NA_integer_
  treedepth_hits <- NA_integer_
  if (!is.null(sampler)) {
    divergences <- sum(vapply(sampler, function(x) sum(x[, "divergent__"] > 0), numeric(1)))
    treedepth_hits <- sum(vapply(sampler, function(x) sum(x[, "treedepth__"] >= max_treedepth), numeric(1)))
  }
  tibble::tibble(
    max_Rhat = max(draw_diag$rhat, na.rm = TRUE),
    min_bulk_ESS = min(draw_diag$ess_bulk, na.rm = TRUE),
    min_tail_ESS = min(draw_diag$ess_tail, na.rm = TRUE),
    divergent_transitions = divergences,
    max_treedepth_hits = treedepth_hits
  )
}

validate_bayes_sens_TE_identity <- function(...) {
  objects <- list(...)
  purrr::imap_dfr(objects, function(obj, name) {
    if (is.null(obj$draws_wide)) return(tibble::tibble(component = name, max_abs_TE_gap = NA_real_, passed = FALSE))
    gap <- obj$draws_wide$TE - obj$draws_wide$direct_effect - obj$draws_wide$total_IIE
    tibble::tibble(component = name, max_abs_TE_gap = max(abs(gap), na.rm = TRUE), passed = max(abs(gap), na.rm = TRUE) < 1e-8)
  })
}

bayes_sens_A_identity_columns <- function() {
  c(
    "direct_effect",
    "total_IIE",
    "TE",
    "PM",
    "alpha_PC1_R",
    "alpha_PC2_R",
    "alpha_PC3",
    "beta_PC1_R",
    "beta_PC2_R",
    "beta_PC3",
    "IIE_PC1_R",
    "IIE_PC2_R",
    "IIE_PC3"
  )
}

validate_bayesian_A_rho0_draw_identity <- function(A_result,
                                                   baseline_quantities,
                                                   tolerance = 1e-12) {
  cols <- bayes_sens_A_identity_columns()
  A0 <- A_result$draws_wide |>
    dplyr::filter(rho_1 == 0, rho_2 == 0, rho_3 == 0) |>
    dplyr::select(
      scenario,
      mediator_varied,
      rho,
      rho_1,
      rho_2,
      rho_3,
      draw,
      exposure,
      contrast,
      dplyr::all_of(cols)
    )
  B0 <- baseline_quantities$wide |>
    dplyr::select(draw, exposure, contrast, dplyr::all_of(cols)) |>
    dplyr::rename_with(~ paste0("baseline_", .x), dplyr::all_of(cols))

  joined <- A0 |>
    dplyr::left_join(B0, by = c("draw", "exposure", "contrast"))

  detail <- purrr::map_dfr(cols, function(col) {
    diff <- joined[[col]] - joined[[paste0("baseline_", col)]]
    tibble::tibble(
      scenario = joined$scenario,
      mediator_varied = joined$mediator_varied,
      draw = joined$draw,
      exposure = joined$exposure,
      contrast = joined$contrast,
      quantity = col,
      A_value = joined[[col]],
      baseline_value = joined[[paste0("baseline_", col)]],
      diff = diff,
      abs_diff = abs(diff)
    )
  })

  max_abs_diff <- max(detail$abs_diff, na.rm = TRUE)
  list(
    passed = is.finite(max_abs_diff) && max_abs_diff <= tolerance,
    tolerance = tolerance,
    max_abs_draw_diff = max_abs_diff,
    max_discrepancy = detail |>
      dplyr::arrange(dplyr::desc(abs_diff)) |>
      dplyr::slice(1),
    detail = detail
  )
}

validate_bayesian_A_rho0_scenarios_identical <- function(A_common,
                                                         A_one_at_a_time,
                                                         tolerance = 1e-12) {
  cols <- bayes_sens_A_identity_columns()
  common0 <- A_common$draws_wide |>
    dplyr::filter(rho_1 == 0, rho_2 == 0, rho_3 == 0) |>
    dplyr::select(draw, exposure, contrast, dplyr::all_of(cols)) |>
    dplyr::rename_with(~ paste0("common_", .x), dplyr::all_of(cols))
  one0 <- A_one_at_a_time$draws_wide |>
    dplyr::filter(rho_1 == 0, rho_2 == 0, rho_3 == 0) |>
    dplyr::select(
      scenario,
      mediator_varied,
      rho,
      rho_1,
      rho_2,
      rho_3,
      draw,
      exposure,
      contrast,
      dplyr::all_of(cols)
    )
  joined <- one0 |>
    dplyr::left_join(common0, by = c("draw", "exposure", "contrast"))

  detail <- purrr::map_dfr(cols, function(col) {
    diff <- joined[[col]] - joined[[paste0("common_", col)]]
    tibble::tibble(
      mediator_varied = joined$mediator_varied,
      draw = joined$draw,
      exposure = joined$exposure,
      contrast = joined$contrast,
      quantity = col,
      one_at_a_time_value = joined[[col]],
      common_value = joined[[paste0("common_", col)]],
      diff = diff,
      abs_diff = abs(diff)
    )
  })
  max_abs_diff <- max(detail$abs_diff, na.rm = TRUE)

  list(
    passed = is.finite(max_abs_diff) && max_abs_diff <= tolerance &&
      identical(A_common$residual_draws$sampled_indices, A_one_at_a_time$residual_draws$sampled_indices),
    tolerance = tolerance,
    sampled_indices_identical = identical(
      A_common$residual_draws$sampled_indices,
      A_one_at_a_time$residual_draws$sampled_indices
    ),
    max_abs_draw_diff = max_abs_diff,
    max_discrepancy = detail |>
      dplyr::arrange(dplyr::desc(abs_diff)) |>
      dplyr::slice(1),
    detail = detail
  )
}

validate_bayesian_A_rho0_against_baseline <- function(A_result,
                                                      baseline_quantities,
                                                      tolerance = 1e-12) {
  A0_summary <- A_result$summary |>
    dplyr::filter(rho_1 == 0, rho_2 == 0, rho_3 == 0) |>
    dplyr::select(
      scenario,
      mediator_varied,
      contrast,
      exposure,
      A_direct_mean = direct_effect_mean,
      A_total_IIE_mean = total_IIE_mean,
      A_TE_mean = TE_mean,
      A_PM_mean = PM_mean
    )
  B_summary <- baseline_quantities$decomposition_summary |>
    dplyr::select(
      contrast,
      exposure,
      baseline_direct_mean = direct_effect_mean,
      baseline_total_IIE_mean = total_IIE_mean,
      baseline_TE_mean = TE_mean,
      baseline_PM_mean = PM_mean
    )
  decomp <- A0_summary |>
    dplyr::left_join(B_summary, by = c("contrast", "exposure")) |>
    dplyr::mutate(
      direct_mean_diff = A_direct_mean - baseline_direct_mean,
      total_IIE_mean_diff = A_total_IIE_mean - baseline_total_IIE_mean,
      TE_mean_diff = A_TE_mean - baseline_TE_mean,
      PM_mean_diff = A_PM_mean - baseline_PM_mean
    )

  A0_path <- A_result$path |>
    dplyr::filter(rho_1 == 0, rho_2 == 0, rho_3 == 0) |>
    dplyr::select(
      scenario,
      mediator_varied,
      contrast,
      exposure,
      mediator,
      A_alpha_mean = alpha_mean,
      A_beta_mean = beta_mean,
      A_IIE_mean = IIE_mean
    )
  B_path <- baseline_quantities$path_summary |>
    dplyr::select(
      contrast,
      exposure,
      mediator,
      baseline_alpha_mean = alpha_mean,
      baseline_beta_mean = beta_mean,
      baseline_IIE_mean = IIE_mean
    )
  path <- A0_path |>
    dplyr::left_join(B_path, by = c("contrast", "exposure", "mediator")) |>
    dplyr::mutate(
      alpha_mean_diff = A_alpha_mean - baseline_alpha_mean,
      beta_mean_diff = A_beta_mean - baseline_beta_mean,
      IIE_mean_diff = A_IIE_mean - baseline_IIE_mean
    )

  max_abs_diff <- max(abs(c(
    decomp$direct_mean_diff,
    decomp$total_IIE_mean_diff,
    decomp$TE_mean_diff,
    decomp$PM_mean_diff,
    path$alpha_mean_diff,
    path$beta_mean_diff,
    path$IIE_mean_diff
  )), na.rm = TRUE)
  draw_identity <- validate_bayesian_A_rho0_draw_identity(
    A_result = A_result,
    baseline_quantities = baseline_quantities,
    tolerance = tolerance
  )

  list(
    passed = isTRUE(draw_identity$passed),
    tolerance = tolerance,
    max_abs_mean_diff = max_abs_diff,
    max_abs_draw_diff = draw_identity$max_abs_draw_diff,
    max_draw_discrepancy = draw_identity$max_discrepancy,
    decomposition = decomp,
    path = path,
    draw_identity = draw_identity
  )
}
