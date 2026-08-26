# Sensitivity B: frequentist latent-binary-U sensitivity analysis.
#
# U is not imputed. For fixed sensitivity parameters, the observed-data
# likelihood is
#   L_i = sum_{u=0}^1 P(U_i = u | X_i)
#         f(M_i | X_i, C_i, U_i = u) f(Y_i | X_i, M_i, C_i, U_i = u).
# lambda_U and delta_U are fixed by standardized user-facing parameters. The
# component means use the sign convention M = ... - lambda_U * U + eps_M and
# Y = ... - delta_U * U + eta. Regression coefficients and residual
# variances/covariances are estimated by ML.

check_sensitivity_B_packages <- function() {
  required <- c("dplyr", "tidyr", "tibble", "purrr", "ggplot2")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Sensitivity B requires installed R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

sensitivity_B_expit <- function(x) {
  stats::plogis(x)
}

sensitivity_B_logit <- function(p) {
  if (!is.finite(p) || p <= 0 || p >= 1) {
    stop("p_U must be strictly between 0 and 1.", call. = FALSE)
  }
  stats::qlogis(p)
}

sensitivity_B_contrast_label <- function(exposure) {
  dplyr::case_when(
    exposure == "X_adjacent" ~ "adjacent vs inside",
    exposure == "X_far" ~ "far vs inside",
    TRUE ~ exposure
  )
}

sensitivity_B_logsumexp2 <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

sensitivity_B_pack_chol <- function(Sigma) {
  L <- t(chol(Sigma))
  c(log(L[1, 1]), L[2, 1], log(L[2, 2]), L[3, 1], L[3, 2], log(L[3, 3]))
}

sensitivity_B_unpack_chol <- function(par) {
  L <- matrix(0, nrow = 3, ncol = 3)
  L[1, 1] <- exp(par[[1]])
  L[2, 1] <- par[[2]]
  L[2, 2] <- exp(par[[3]])
  L[3, 1] <- par[[4]]
  L[3, 2] <- par[[5]]
  L[3, 3] <- exp(par[[6]])
  L %*% t(L)
}

sensitivity_B_logdmvnorm3 <- function(residual, Sigma) {
  if (any(!is.finite(Sigma))) {
    return(rep(-Inf, nrow(residual)))
  }
  chol_u <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(chol_u)) {
    return(rep(-Inf, nrow(residual)))
  }
  solved <- backsolve(chol_u, t(residual), transpose = TRUE)
  quad <- colSums(solved^2)
  logdet <- 2 * sum(log(diag(chol_u)))
  if (!is.finite(logdet) || any(!is.finite(quad))) {
    return(rep(-Inf, nrow(residual)))
  }
  -0.5 * (3 * log(2 * pi) + logdet + quad)
}

sensitivity_B_baseline_fit <- function(data,
                                       exposures = c("X_adjacent", "X_far"),
                                       outcome = "Y",
                                       mediators = c("PC1_R", "PC2_R", "PC3"),
                                       covariates = c("x_coord", "y_coord")) {
  X_m <- stats::model.matrix(
    stats::as.formula(paste("~", paste(c(exposures, covariates), collapse = " + "))),
    data = data
  )
  M <- as.matrix(data[, mediators, drop = FALSE])
  B_m <- solve(crossprod(X_m), crossprod(X_m, M))
  E_m <- M - X_m %*% B_m
  Sigma_M <- crossprod(E_m) / nrow(E_m)

  X_y <- stats::model.matrix(
    stats::as.formula(paste(outcome, "~", paste(c(exposures, mediators, covariates), collapse = " + "))),
    data = data
  )
  y <- data[[outcome]]
  b_y <- solve(crossprod(X_y), crossprod(X_y, y))
  e_y <- as.numeric(y - X_y %*% b_y)
  sigma_y <- sqrt(mean(e_y^2))

  X_t <- stats::model.matrix(
    stats::as.formula(paste(outcome, "~", paste(c(exposures, covariates), collapse = " + "))),
    data = data
  )
  b_t <- solve(crossprod(X_t), crossprod(X_t, y))
  e_t <- as.numeric(y - X_t %*% b_t)
  sigma_t <- sqrt(mean(e_t^2))

  decomp <- decompose_linear_multix(
    data = data,
    exposures = exposures,
    outcome = outcome,
    mediators = mediators,
    covariates = covariates
  )

  list(
    X_m = X_m,
    M = M,
    X_y = X_y,
    y = y,
    X_t = X_t,
    B_m = B_m,
    Sigma_M = Sigma_M,
    b_y = as.numeric(b_y),
    sigma_y = sigma_y,
    b_t = as.numeric(b_t),
    sigma_t = sigma_t,
    mediator_terms = colnames(X_m),
    outcome_terms = colnames(X_y),
    total_terms = colnames(X_t),
    decomposition = decomp,
    summary = decomp$summary |> dplyr::mutate(contrast = sensitivity_B_contrast_label(exposure)),
    path = decomp$path |> dplyr::mutate(contrast = sensitivity_B_contrast_label(exposure))
  )
}

sensitivity_B_make_fixed_params <- function(data,
                                            p_U = 0.3,
                                            OR_U_adj = 1,
                                            OR_U_far = 1,
                                            s_M = c(0, 0, 0),
                                            s_Y = 0,
                                            mediators = c("PC1_R", "PC2_R", "PC3"),
                                            outcome = "Y") {
  if (length(s_M) == 1) {
    s_M <- rep(s_M, length(mediators))
  }
  if (length(s_M) != length(mediators)) {
    stop("s_M must be length 1 or length equal to number of mediators.", call. = FALSE)
  }
  if (OR_U_adj <= 0 || OR_U_far <= 0) {
    stop("OR_U_adj and OR_U_far must be positive.", call. = FALSE)
  }

  sd_M <- vapply(mediators, function(m) stats::sd(data[[m]], na.rm = TRUE), numeric(1))
  sd_Y <- stats::sd(data[[outcome]], na.rm = TRUE)

  list(
    p_U = p_U,
    gamma0 = sensitivity_B_logit(p_U),
    gamma_adj = log(OR_U_adj),
    gamma_far = log(OR_U_far),
    OR_U_adj = OR_U_adj,
    OR_U_far = OR_U_far,
    s_M = stats::setNames(as.numeric(s_M), mediators),
    s_Y = s_Y,
    lambda_U = stats::setNames(as.numeric(s_M) * sd_M, mediators),
    delta_U = s_Y * sd_Y,
    sd_M = sd_M,
    sd_Y = sd_Y
  )
}

sensitivity_B_pi <- function(data, fixed) {
  sensitivity_B_expit(fixed$gamma0 + fixed$gamma_adj * data$X_adjacent + fixed$gamma_far * data$X_far)
}

sensitivity_B_is_baseline_equivalent <- function(fixed, tolerance = 0) {
  no_x_u <- abs(fixed$OR_U_adj - 1) <= tolerance && abs(fixed$OR_U_far - 1) <= tolerance
  no_m_u <- all(abs(fixed$lambda_U) <= tolerance)
  no_y_u <- abs(fixed$delta_U) <= tolerance

  (no_m_u && no_y_u) || (no_x_u && (no_m_u || no_y_u))
}

sensitivity_B_pack_joint <- function(B_m, Sigma_M, b_y, sigma_y) {
  c(as.vector(B_m), sensitivity_B_pack_chol(Sigma_M), b_y, log(sigma_y))
}

sensitivity_B_positive_definite_start <- function(Sigma, min_eigen = 1e-6) {
  Sigma <- (Sigma + t(Sigma)) / 2
  eig <- eigen(Sigma, symmetric = TRUE)
  values <- pmax(eig$values, min_eigen)
  Sigma_pd <- eig$vectors %*% diag(values, nrow = length(values)) %*% t(eig$vectors)
  (Sigma_pd + t(Sigma_pd)) / 2
}

sensitivity_B_projected_u_variance <- function(X, pi_u) {
  pi_coef <- as.numeric(solve(crossprod(X), crossprod(X, pi_u)))
  pi_hat <- as.numeric(X %*% pi_coef)
  list(
    coef = pi_coef,
    variance = mean(pi_u * (1 - pi_u) + (pi_u - pi_hat)^2)
  )
}

sensitivity_B_initial_values <- function(data, baseline, fixed) {
  pi_u <- sensitivity_B_pi(data, fixed)

  u_m <- sensitivity_B_projected_u_variance(baseline$X_m, pi_u)
  B_m_start <- baseline$B_m
  for (j in seq_len(ncol(B_m_start))) {
    B_m_start[, j] <- B_m_start[, j] + u_m$coef * fixed$lambda_U[[j]]
  }

  u_y <- sensitivity_B_projected_u_variance(baseline$X_y, pi_u)
  b_y_start <- baseline$b_y + u_y$coef * fixed$delta_U

  lambda <- as.numeric(fixed$lambda_U)
  beta_terms <- match(colnames(baseline$M), baseline$outcome_terms)
  if (any(abs(lambda) > 0) && abs(fixed$delta_U) > 0) {
    # Under the requested sign convention, latent U contributes approximately
    # Var(U residual) * lambda_U * delta_U to the observed M-Y residual
    # covariance. Start beta on the structural scale rather than at the
    # baseline observed regression coefficient.
    beta_shift <- as.numeric(qr.solve(baseline$Sigma_M, lambda * fixed$delta_U * u_m$variance))
    b_y_start[beta_terms] <- b_y_start[beta_terms] - beta_shift
  }

  u_t <- sensitivity_B_projected_u_variance(baseline$X_t, pi_u)
  b_t_start <- baseline$b_t + u_t$coef * fixed$delta_U

  Sigma_M_start <- baseline$Sigma_M - u_m$variance * tcrossprod(lambda)
  Sigma_M_start <- sensitivity_B_positive_definite_start(
    Sigma_M_start,
    min_eigen = max(1e-6, min(diag(baseline$Sigma_M)) * 1e-8)
  )

  sigma_y2_start <- baseline$sigma_y^2 - u_y$variance * fixed$delta_U^2
  sigma_y_start <- sqrt(max(sigma_y2_start, baseline$sigma_y^2 * 1e-4, 1e-8))

  sigma_t2_start <- baseline$sigma_t^2 - u_t$variance * fixed$delta_U^2
  sigma_t_start <- sqrt(max(sigma_t2_start, baseline$sigma_t^2 * 1e-4, 1e-8))

  list(
    B_m = B_m_start,
    b_y = b_y_start,
    b_t = b_t_start,
    Sigma_M = Sigma_M_start,
    sigma_y = sigma_y_start,
    sigma_t = sigma_t_start,
    projected_u_variance = c(
      mediator = u_m$variance,
      outcome = u_y$variance,
      total = u_t$variance
    )
  )
}

sensitivity_B_unpack_joint <- function(par, p_m, k_m, p_y) {
  n_B <- p_m * k_m
  B_m <- matrix(par[seq_len(n_B)], nrow = p_m, ncol = k_m)
  Sigma_M <- sensitivity_B_unpack_chol(par[n_B + seq_len(6)])
  y_start <- n_B + 6
  b_y <- par[y_start + seq_len(p_y)]
  sigma_y <- exp(par[y_start + p_y + 1])
  list(B_m = B_m, Sigma_M = Sigma_M, b_y = b_y, sigma_y = sigma_y)
}

sensitivity_B_negloglik_joint <- function(par, data, baseline, fixed) {
  unpacked <- sensitivity_B_unpack_joint(
    par = par,
    p_m = ncol(baseline$X_m),
    k_m = ncol(baseline$M),
    p_y = ncol(baseline$X_y)
  )
  if (!is.finite(unpacked$sigma_y) || unpacked$sigma_y <= 0 || any(!is.finite(unpacked$B_m)) || any(!is.finite(unpacked$b_y))) {
    return(.Machine$double.xmax^0.25)
  }

  pi_u <- sensitivity_B_pi(data, fixed)
  mu_m0 <- baseline$X_m %*% unpacked$B_m
  mu_m1 <- sweep(mu_m0, 2, fixed$lambda_U, "-")
  res_m0 <- baseline$M - mu_m0
  res_m1 <- baseline$M - mu_m1

  mu_y0 <- as.numeric(baseline$X_y %*% unpacked$b_y)
  mu_y1 <- mu_y0 - fixed$delta_U
  y <- baseline$y

  log_m0 <- sensitivity_B_logdmvnorm3(res_m0, unpacked$Sigma_M)
  log_m1 <- sensitivity_B_logdmvnorm3(res_m1, unpacked$Sigma_M)
  log_y0 <- stats::dnorm(y, mean = mu_y0, sd = unpacked$sigma_y, log = TRUE)
  log_y1 <- stats::dnorm(y, mean = mu_y1, sd = unpacked$sigma_y, log = TRUE)

  ll0 <- log1p(-pi_u) + log_m0 + log_y0
  ll1 <- log(pi_u) + log_m1 + log_y1
  nll <- -sum(sensitivity_B_logsumexp2(ll0, ll1))
  if (!is.finite(nll)) .Machine$double.xmax^0.25 else nll
}

sensitivity_B_pack_total <- function(b_t, sigma_t) {
  c(b_t, log(sigma_t))
}

sensitivity_B_unpack_total <- function(par) {
  list(b_t = par[-length(par)], sigma_t = exp(par[[length(par)]]))
}

sensitivity_B_negloglik_total <- function(par, data, baseline, fixed, u_effect_y = fixed$delta_U) {
  unpacked <- sensitivity_B_unpack_total(par)
  if (!is.finite(unpacked$sigma_t) || unpacked$sigma_t <= 0 || any(!is.finite(unpacked$b_t))) {
    return(.Machine$double.xmax^0.25)
  }
  pi_u <- sensitivity_B_pi(data, fixed)
  mu0 <- as.numeric(baseline$X_t %*% unpacked$b_t)
  mu1 <- mu0 - u_effect_y
  y <- baseline$y
  ll0 <- log1p(-pi_u) + stats::dnorm(y, mean = mu0, sd = unpacked$sigma_t, log = TRUE)
  ll1 <- log(pi_u) + stats::dnorm(y, mean = mu1, sd = unpacked$sigma_t, log = TRUE)
  nll <- -sum(sensitivity_B_logsumexp2(ll0, ll1))
  if (!is.finite(nll)) .Machine$double.xmax^0.25 else nll
}

sensitivity_B_build_result <- function(baseline,
                                       fixed,
                                       scenario,
                                       B_m,
                                       Sigma_M,
                                       b_y,
                                       sigma_y,
                                       b_t,
                                       sigma_t,
                                       joint_convergence = 0L,
                                       joint_value = NA_real_,
                                       total_convergence = 0L,
                                       total_value = NA_real_,
                                       joint_elapsed_seconds = NA_real_,
                                       total_elapsed_seconds = NA_real_,
                                       used_baseline_equivalence = FALSE,
                                       total_u_effect_y = fixed$delta_U,
                                       raw_optim = NULL) {
  rownames(B_m) <- baseline$mediator_terms
  colnames(B_m) <- colnames(baseline$M)
  names(b_y) <- baseline$outcome_terms
  names(b_t) <- baseline$total_terms

  alpha <- B_m[c("X_adjacent", "X_far"), , drop = FALSE]
  beta <- b_y[colnames(baseline$M)]
  direct <- b_y[c("X_adjacent", "X_far")]
  TE_reduced_form <- b_t[c("X_adjacent", "X_far")]
  IIE <- sweep(alpha, 2, beta, "*")
  total_IIE <- rowSums(IIE)
  TE_decomp <- direct + total_IIE

  path <- as.data.frame(alpha) |>
    tibble::rownames_to_column("exposure") |>
    tidyr::pivot_longer(cols = dplyr::all_of(colnames(baseline$M)), names_to = "mediator", values_to = "alpha") |>
    dplyr::left_join(tibble::tibble(mediator = names(beta), beta = as.numeric(beta)), by = "mediator") |>
    dplyr::mutate(
      scenario = scenario,
      contrast = sensitivity_B_contrast_label(exposure),
      IIE = alpha * beta,
      p_U = fixed$p_U,
      OR_U_adj = fixed$OR_U_adj,
      OR_U_far = fixed$OR_U_far,
      s_Y = fixed$s_Y,
      s_M1 = fixed$s_M[[1]],
      s_M2 = fixed$s_M[[2]],
      s_M3 = fixed$s_M[[3]],
      lambda_U = fixed$lambda_U[mediator],
      delta_U = fixed$delta_U
    ) |>
    dplyr::select(scenario, p_U, OR_U_adj, OR_U_far, s_M1, s_M2, s_M3, s_Y, contrast, exposure, mediator, alpha, beta, IIE, lambda_U, delta_U)

  exposure_vec <- c("X_adjacent", "X_far")
  summary <- tibble::tibble(
    scenario = scenario,
    p_U = fixed$p_U,
    OR_U_adj = fixed$OR_U_adj,
    OR_U_far = fixed$OR_U_far,
    s_M1 = fixed$s_M[[1]],
    s_M2 = fixed$s_M[[2]],
    s_M3 = fixed$s_M[[3]],
    s_Y = fixed$s_Y,
    exposure = exposure_vec,
    contrast = sensitivity_B_contrast_label(exposure_vec),
    TE = TE_decomp[exposure_vec],
    TE_decomp = TE_decomp[exposure_vec],
    TE_reduced_form = TE_reduced_form[exposure_vec],
    direct_effect = direct[exposure_vec],
    total_IIE = total_IIE[exposure_vec],
    PM = total_IIE[exposure_vec] / TE_decomp[exposure_vec],
    delta_U = fixed$delta_U,
    total_u_effect_y = total_u_effect_y
  )

  convergence <- tibble::tibble(
    scenario = scenario,
    p_U = fixed$p_U,
    OR_U_adj = fixed$OR_U_adj,
    OR_U_far = fixed$OR_U_far,
    s_M1 = fixed$s_M[[1]],
    s_M2 = fixed$s_M[[2]],
    s_M3 = fixed$s_M[[3]],
    s_Y = fixed$s_Y,
    joint_convergence = joint_convergence,
    joint_value = joint_value,
    total_convergence = total_convergence,
    total_value = total_value,
    joint_elapsed_seconds = joint_elapsed_seconds,
    total_elapsed_seconds = total_elapsed_seconds,
    used_baseline_equivalence = used_baseline_equivalence,
    total_u_effect_y = total_u_effect_y
  )

  list(
    fixed = fixed,
    estimates = list(
      joint = list(B_m = B_m, Sigma_M = Sigma_M, b_y = b_y, sigma_y = sigma_y),
      total = list(b_t = b_t, sigma_t = sigma_t)
    ),
    path = path,
    summary = summary,
    convergence = convergence,
    raw_optim = raw_optim
  )
}

fit_sensitivity_B_one <- function(data,
                                  baseline,
                                  scenario,
                                  p_U = 0.3,
                                  OR_U_adj = 1,
                                  OR_U_far = 1,
                                  s_M = c(0, 0, 0),
                                  s_Y = 0,
                                  maxit = 1000,
                                  reltol = 1e-8) {
  fixed <- sensitivity_B_make_fixed_params(
    data = data,
    p_U = p_U,
    OR_U_adj = OR_U_adj,
    OR_U_far = OR_U_far,
    s_M = s_M,
    s_Y = s_Y
  )

  if (sensitivity_B_is_baseline_equivalent(fixed)) {
    return(sensitivity_B_build_result(
      baseline = baseline,
      fixed = fixed,
      scenario = scenario,
      B_m = baseline$B_m,
      Sigma_M = baseline$Sigma_M,
      b_y = baseline$b_y,
      sigma_y = baseline$sigma_y,
      b_t = baseline$b_t,
      sigma_t = baseline$sigma_t,
      used_baseline_equivalence = TRUE,
      total_u_effect_y = fixed$delta_U + sum(baseline$b_y[match(colnames(baseline$M), baseline$outcome_terms)] * as.numeric(fixed$lambda_U))
    ))
  }

  starts <- sensitivity_B_initial_values(data, baseline, fixed)
  joint_start <- sensitivity_B_pack_joint(
    B_m = starts$B_m,
    Sigma_M = starts$Sigma_M,
    b_y = starts$b_y,
    sigma_y = starts$sigma_y
  )

  joint_time <- system.time({
    fit_joint <- stats::optim(
      par = joint_start,
      fn = sensitivity_B_negloglik_joint,
      data = data,
      baseline = baseline,
      fixed = fixed,
      method = "BFGS",
      control = list(maxit = maxit, reltol = reltol)
    )
  })

  joint <- sensitivity_B_unpack_joint(
    fit_joint$par,
    p_m = ncol(baseline$X_m),
    k_m = ncol(baseline$M),
    p_y = ncol(baseline$X_y)
  )
  beta_for_total <- joint$b_y[match(colnames(baseline$M), baseline$outcome_terms)]
  total_u_effect_y <- fixed$delta_U + sum(beta_for_total * as.numeric(fixed$lambda_U))

  pi_u <- sensitivity_B_pi(data, fixed)
  u_t <- sensitivity_B_projected_u_variance(baseline$X_t, pi_u)
  b_t_start <- baseline$b_t + u_t$coef * total_u_effect_y
  sigma_t2_start <- baseline$sigma_t^2 - u_t$variance * total_u_effect_y^2
  sigma_t_start <- sqrt(max(sigma_t2_start, baseline$sigma_t^2 * 1e-4, 1e-8))
  total_start <- sensitivity_B_pack_total(b_t_start, sigma_t_start)

  total_time <- system.time({
    fit_total <- stats::optim(
      par = total_start,
      fn = sensitivity_B_negloglik_total,
      data = data,
      baseline = baseline,
      fixed = fixed,
      u_effect_y = total_u_effect_y,
      method = "BFGS",
      control = list(maxit = maxit, reltol = reltol)
    )
  })

  total <- sensitivity_B_unpack_total(fit_total$par)

  sensitivity_B_build_result(
    baseline = baseline,
    fixed = fixed,
    scenario = scenario,
    B_m = joint$B_m,
    Sigma_M = joint$Sigma_M,
    b_y = joint$b_y,
    sigma_y = joint$sigma_y,
    b_t = total$b_t,
    sigma_t = total$sigma_t,
    joint_convergence = fit_joint$convergence,
    joint_value = fit_joint$value,
    total_convergence = fit_total$convergence,
    total_value = fit_total$value,
    joint_elapsed_seconds = unname(joint_time[["elapsed"]]),
    total_elapsed_seconds = unname(total_time[["elapsed"]]),
    used_baseline_equivalence = FALSE,
    total_u_effect_y = total_u_effect_y,
    raw_optim = list(joint = fit_joint, total = fit_total)
  )
}

sensitivity_B_add_changes <- function(summary, path, baseline) {
  base_summary <- baseline$summary |>
    dplyr::select(exposure, baseline_TE = TE, baseline_direct = NDE, baseline_total_IIE = NIE, baseline_PM = PM)
  base_path <- baseline$path |>
    dplyr::select(exposure, mediator, baseline_alpha = alpha_X_to_M, baseline_beta = beta_M_to_Y, baseline_IIE = indirect_component)

  list(
    summary = summary |>
      dplyr::left_join(base_summary, by = "exposure") |>
      dplyr::mutate(
        TE_change = TE - baseline_TE,
        TE_reduced_form_change = TE_reduced_form - baseline_TE,
        direct_change = direct_effect - baseline_direct,
        total_IIE_change = total_IIE - baseline_total_IIE,
        PM_change = PM - baseline_PM,
        total_IIE_relative_change = total_IIE_change / baseline_total_IIE
      ),
    path = path |>
      dplyr::left_join(base_path, by = c("exposure", "mediator")) |>
      dplyr::mutate(
        alpha_change = alpha - baseline_alpha,
        beta_change = beta - baseline_beta,
        IIE_change = IIE - baseline_IIE,
        IIE_relative_change = IIE_change / baseline_IIE
      )
  )
}

make_sensitivity_B_grid <- function(scenario,
                                    p_U_grid = 0.3,
                                    OR_UX_grid = c(1, 1.5, 2, 3),
                                    s_M_grid = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5),
                                    s_Y_grid = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5)) {
  scenario <- match.arg(scenario, c("B1_XY", "B2_XM", "B3_MY", "B4_global"))

  if (scenario == "B1_XY") {
    return(expand.grid(p_U = p_U_grid, OR_UX = OR_UX_grid, s_Y = s_Y_grid) |>
      tibble::as_tibble() |>
      dplyr::mutate(s_M = 0))
  }
  if (scenario == "B2_XM") {
    return(expand.grid(p_U = p_U_grid, OR_UX = OR_UX_grid, s_M = s_M_grid) |>
      tibble::as_tibble() |>
      dplyr::mutate(s_Y = 0))
  }
  if (scenario == "B3_MY") {
    return(expand.grid(p_U = p_U_grid, s_M = s_M_grid, s_Y = s_Y_grid) |>
      tibble::as_tibble() |>
      dplyr::mutate(OR_UX = 1))
  }

  expand.grid(p_U = p_U_grid, OR_UX = OR_UX_grid, s_M = s_M_grid, s_Y = s_Y_grid) |>
    tibble::as_tibble()
}

run_sensitivity_B_grid <- function(data,
                                   baseline,
                                   scenario,
                                   grid,
                                   maxit = 1000,
                                   reltol = 1e-8) {
  results <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    row <- grid[i, ]
    results[[i]] <- tryCatch(
      fit_sensitivity_B_one(
        data = data,
        baseline = baseline,
        scenario = scenario,
        p_U = row$p_U,
        OR_U_adj = row$OR_UX,
        OR_U_far = row$OR_UX,
        s_M = row$s_M,
        s_Y = row$s_Y,
        maxit = maxit,
        reltol = reltol
      ),
      error = function(e) {
        list(
          path = tibble::tibble(),
          summary = tibble::tibble(),
          convergence = dplyr::bind_cols(row, tibble::tibble(
            scenario = scenario,
            joint_convergence = NA_integer_,
            joint_value = NA_real_,
            total_convergence = NA_integer_,
            total_value = NA_real_,
            joint_elapsed_seconds = NA_real_,
            total_elapsed_seconds = NA_real_,
            failure_reason = conditionMessage(e)
          ))
        )
      }
    )
  }

  path <- dplyr::bind_rows(purrr::map(results, "path"))
  summary <- dplyr::bind_rows(purrr::map(results, "summary"))
  convergence <- dplyr::bind_rows(purrr::map(results, "convergence"))
  changed <- sensitivity_B_add_changes(summary, path, baseline)

  list(
    grid = grid,
    raw = results,
    summary = changed$summary,
    path = changed$path,
    convergence = convergence
  )
}

validate_sensitivity_B_null <- function(data,
                                        baseline,
                                        tolerance = 1e-5,
                                        maxit = 1000,
                                        reltol = 1e-9) {
  null_fit <- fit_sensitivity_B_one(
    data = data,
    baseline = baseline,
    scenario = "null_all",
    p_U = 0.3,
    OR_U_adj = 1,
    OR_U_far = 1,
    s_M = c(0, 0, 0),
    s_Y = 0,
    maxit = maxit,
    reltol = reltol
  )
  changed <- sensitivity_B_add_changes(null_fit$summary, null_fit$path, baseline)

  list(
    summary = changed$summary,
    path = changed$path,
    convergence = null_fit$convergence,
    passed = max(abs(changed$summary$total_IIE_change), na.rm = TRUE) <= tolerance &&
      max(abs(changed$summary$direct_change), na.rm = TRUE) <= tolerance &&
      max(abs(changed$summary$TE_change), na.rm = TRUE) <= tolerance &&
      max(abs(changed$path$IIE_change), na.rm = TRUE) <= tolerance,
    tolerance = tolerance
  )
}

validate_sensitivity_B_scenario_nulls <- function(data,
                                                  baseline,
                                                  tolerance = 1e-4,
                                                  p_U = 0.3,
                                                  OR_test = 2,
                                                  s_M_test = 0.3,
                                                  s_Y_test = 0.3,
                                                  maxit = 1000,
                                                  reltol = 1e-8) {
  checks <- tibble::tribble(
    ~check, ~scenario, ~OR_U_adj, ~OR_U_far, ~s_M, ~s_Y,
    "B1_OR_UX_eq_1", "B1_XY", 1, 1, 0, s_Y_test,
    "B1_s_Y_eq_0", "B1_XY", OR_test, OR_test, 0, 0,
    "B2_OR_UX_eq_1", "B2_XM", 1, 1, s_M_test, 0,
    "B2_s_M_eq_0", "B2_XM", OR_test, OR_test, 0, 0,
    "B3_s_M_eq_0", "B3_MY", 1, 1, 0, s_Y_test,
    "B3_s_Y_eq_0", "B3_MY", 1, 1, s_M_test, 0,
    "B4_all_null", "B4_global", 1, 1, 0, 0
  )

  purrr::pmap_dfr(checks, function(check, scenario, OR_U_adj, OR_U_far, s_M, s_Y) {
    fit <- tryCatch(
      fit_sensitivity_B_one(
        data = data,
        baseline = baseline,
        scenario = scenario,
        p_U = p_U,
        OR_U_adj = OR_U_adj,
        OR_U_far = OR_U_far,
        s_M = s_M,
        s_Y = s_Y,
        maxit = maxit,
        reltol = reltol
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(tibble::tibble(
        check = check,
        scenario = scenario,
        passed = FALSE,
        max_abs_summary_change = NA_real_,
        max_abs_path_change = NA_real_,
        failure_reason = conditionMessage(fit)
      ))
    }

    changed <- sensitivity_B_add_changes(fit$summary, fit$path, baseline)
    summary_change <- max(abs(c(
      changed$summary$TE_change,
      changed$summary$direct_change,
      changed$summary$total_IIE_change,
      changed$summary$PM_change
    )), na.rm = TRUE)
    path_change <- max(abs(c(
      changed$path$alpha_change,
      changed$path$beta_change,
      changed$path$IIE_change
    )), na.rm = TRUE)

    tibble::tibble(
      check = check,
      scenario = scenario,
      passed = summary_change <= tolerance && path_change <= tolerance,
      max_abs_summary_change = summary_change,
      max_abs_path_change = path_change,
      failure_reason = NA_character_
    )
  })
}

sensitivity_B_changed_result <- function(fit, baseline) {
  changed <- sensitivity_B_add_changes(fit$summary, fit$path, baseline)
  list(
    summary = changed$summary,
    path = changed$path,
    convergence = fit$convergence
  )
}

sensitivity_B_compare_changed <- function(left, right) {
  path_joined <- left$path |>
    dplyr::select(exposure, mediator, alpha, beta, IIE) |>
    dplyr::rename(
      alpha_left = alpha,
      beta_left = beta,
      IIE_left = IIE
    ) |>
    dplyr::inner_join(
      right$path |>
        dplyr::select(exposure, mediator, alpha, beta, IIE) |>
        dplyr::rename(
          alpha_right = alpha,
          beta_right = beta,
          IIE_right = IIE
        ),
      by = c("exposure", "mediator")
    )

  summary_joined <- left$summary |>
    dplyr::select(exposure, TE, direct_effect, total_IIE, PM) |>
    dplyr::rename(
      TE_left = TE,
      direct_left = direct_effect,
      total_IIE_left = total_IIE,
      PM_left = PM
    ) |>
    dplyr::inner_join(
      right$summary |>
        dplyr::select(exposure, TE, direct_effect, total_IIE, PM) |>
        dplyr::rename(
          TE_right = TE,
          direct_right = direct_effect,
          total_IIE_right = total_IIE,
          PM_right = PM
        ),
      by = "exposure"
    )

  max(abs(c(
    path_joined$alpha_left - path_joined$alpha_right,
    path_joined$beta_left - path_joined$beta_right,
    path_joined$IIE_left - path_joined$IIE_right,
    summary_joined$TE_left - summary_joined$TE_right,
    summary_joined$direct_left - summary_joined$direct_right,
    summary_joined$total_IIE_left - summary_joined$total_IIE_right,
    summary_joined$PM_left - summary_joined$PM_right
  )), na.rm = TRUE)
}

validate_sensitivity_B_decomposition_identity <- function(result_list, tolerance = 1e-8) {
  purrr::imap_dfr(result_list, function(res, scenario) {
    if (is.null(res) || is.null(res$summary) || nrow(res$summary) == 0) {
      return(tibble::tibble(
        scenario = scenario,
        passed = FALSE,
        max_abs_decomposition_gap = NA_real_,
        failure_reason = "No summary rows."
      ))
    }
    gap <- res$summary$TE - res$summary$direct_effect - res$summary$total_IIE
    tibble::tibble(
      scenario = scenario,
      passed = max(abs(gap), na.rm = TRUE) <= tolerance,
      max_abs_decomposition_gap = max(abs(gap), na.rm = TRUE),
      failure_reason = NA_character_
    )
  })
}

validate_sensitivity_B_path_behavior <- function(result_list,
                                                 tolerance = 1e-4,
                                                 nonzero_tolerance = 1e-4) {
  get_res <- function(name) {
    result_list[[name]]
  }
  path_check <- function(df, check, condition, quantity, expect_changed = FALSE) {
    rows <- df |> dplyr::filter({{ condition }})
    if (nrow(rows) == 0) {
      return(tibble::tibble(
        check = check,
        passed = FALSE,
        max_abs_change = NA_real_,
        failure_reason = "No matching rows."
      ))
    }
    max_change <- max(abs(rows[[quantity]]), na.rm = TRUE)
    tibble::tibble(
      check = check,
      passed = if (expect_changed) max_change > nonzero_tolerance else max_change <= tolerance,
      max_abs_change = max_change,
      failure_reason = NA_character_
    )
  }

  out <- list()
  if (!is.null(get_res("B1_XY"))) {
    b1 <- get_res("B1_XY")$path
    out[["B1_alpha_baseline"]] <- path_check(b1, "B1_XY_alpha_approximately_baseline", TRUE, "alpha_change")
    out[["B1_beta_baseline"]] <- path_check(b1, "B1_XY_beta_approximately_baseline", TRUE, "beta_change")
  }
  if (!is.null(get_res("B2_XM"))) {
    b2 <- get_res("B2_XM")$path
    out[["B2_alpha_changes"]] <- path_check(
      b2,
      "B2_XM_alpha_changes_when_OR_UX_and_s_M_nonzero",
      OR_U_adj != 1 & (s_M1 != 0 | s_M2 != 0 | s_M3 != 0),
      "alpha_change",
      expect_changed = TRUE
    )
    out[["B2_beta_baseline"]] <- path_check(b2, "B2_XM_beta_approximately_baseline", TRUE, "beta_change")
  }
  if (!is.null(get_res("B3_MY"))) {
    b3 <- get_res("B3_MY")$path
    out[["B3_alpha_baseline"]] <- path_check(b3, "B3_MY_alpha_approximately_baseline", TRUE, "alpha_change")
    out[["B3_beta_changes"]] <- path_check(
      b3,
      "B3_MY_beta_changes_when_s_M_and_s_Y_nonzero",
      (s_M1 != 0 | s_M2 != 0 | s_M3 != 0) & s_Y != 0,
      "beta_change",
      expect_changed = TRUE
    )
  }
  if (!is.null(get_res("B4_global"))) {
    b4 <- get_res("B4_global")$path
    out[["B4_alpha_may_change"]] <- path_check(
      b4,
      "B4_global_alpha_can_change",
      OR_U_adj != 1 & (s_M1 != 0 | s_M2 != 0 | s_M3 != 0),
      "alpha_change",
      expect_changed = TRUE
    )
    out[["B4_beta_may_change"]] <- path_check(
      b4,
      "B4_global_beta_can_change",
      (s_M1 != 0 | s_M2 != 0 | s_M3 != 0) & s_Y != 0,
      "beta_change",
      expect_changed = TRUE
    )
  }

  dplyr::bind_rows(out)
}

validate_sensitivity_B_nesting <- function(data,
                                           baseline,
                                           tolerance = 1e-4,
                                           nonzero_tolerance = 1e-4,
                                           p_U = 0.3,
                                           OR_test = 2,
                                           s_M_test = 0.5,
                                           s_Y_test = 0.5,
                                           maxit = 1000,
                                           reltol = 1e-8) {
  fit_changed <- function(...) {
    fit <- fit_sensitivity_B_one(
      data = data,
      baseline = baseline,
      p_U = p_U,
      maxit = maxit,
      reltol = reltol,
      ...
    )
    sensitivity_B_changed_result(fit, baseline)
  }

  baseline_checks <- tibble::tribble(
    ~check, ~scenario, ~OR_U_adj, ~OR_U_far, ~s_M, ~s_Y,
    "B3_OR1_sM_nonzero_sY_nonzero_beta_differs", "B3_MY", 1, 1, s_M_test, s_Y_test,
    "B3_OR1_sM_zero_sY_nonzero_matches_baseline", "B3_MY", 1, 1, 0, s_Y_test,
    "B3_OR1_sM_nonzero_sY_zero_matches_baseline", "B3_MY", 1, 1, s_M_test, 0
  )

  baseline_result <- purrr::pmap_dfr(
    baseline_checks,
    function(check, scenario, OR_U_adj, OR_U_far, s_M, s_Y) {
      res <- tryCatch(
        fit_changed(
          scenario = scenario,
          OR_U_adj = OR_U_adj,
          OR_U_far = OR_U_far,
          s_M = s_M,
          s_Y = s_Y
        ),
        error = function(e) e
      )
      if (inherits(res, "error")) {
        return(tibble::tibble(
          check = check,
          comparison = "baseline",
          passed = FALSE,
          max_abs_beta_change = NA_real_,
          max_abs_IIE_change = NA_real_,
          max_abs_total_IIE_change = NA_real_,
          max_abs_comparison_difference = NA_real_,
          failure_reason = conditionMessage(res)
        ))
      }

      max_beta <- max(abs(res$path$beta_change), na.rm = TRUE)
      max_iie <- max(abs(res$path$IIE_change), na.rm = TRUE)
      max_total_iie <- max(abs(res$summary$total_IIE_change), na.rm = TRUE)
      should_change <- grepl("beta_differs", check, fixed = TRUE)
      passed <- if (should_change) {
        max_beta > nonzero_tolerance && max_iie > nonzero_tolerance
      } else {
        max(abs(c(max_beta, max_iie, max_total_iie)), na.rm = TRUE) <= tolerance
      }

      tibble::tibble(
        check = check,
        comparison = "baseline",
        passed = passed,
        max_abs_beta_change = max_beta,
        max_abs_IIE_change = max_iie,
        max_abs_total_IIE_change = max_total_iie,
        max_abs_comparison_difference = NA_real_,
        failure_reason = NA_character_
      )
    }
  )

  pair_specs <- list(
    list(
      check = "B4_OR1_sM_nonzero_sY_nonzero_matches_B3",
      left = list(scenario = "B4_global", OR_U_adj = 1, OR_U_far = 1, s_M = s_M_test, s_Y = s_Y_test),
      right = list(scenario = "B3_MY", OR_U_adj = 1, OR_U_far = 1, s_M = s_M_test, s_Y = s_Y_test)
    ),
    list(
      check = "B4_ORgt1_sM_zero_sY_nonzero_matches_B1",
      left = list(scenario = "B4_global", OR_U_adj = OR_test, OR_U_far = OR_test, s_M = 0, s_Y = s_Y_test),
      right = list(scenario = "B1_XY", OR_U_adj = OR_test, OR_U_far = OR_test, s_M = 0, s_Y = s_Y_test)
    ),
    list(
      check = "B4_ORgt1_sM_nonzero_sY_zero_matches_B2",
      left = list(scenario = "B4_global", OR_U_adj = OR_test, OR_U_far = OR_test, s_M = s_M_test, s_Y = 0),
      right = list(scenario = "B2_XM", OR_U_adj = OR_test, OR_U_far = OR_test, s_M = s_M_test, s_Y = 0)
    )
  )

  pair_result <- purrr::map_dfr(pair_specs, function(spec) {
    left <- tryCatch(do.call(fit_changed, spec$left), error = function(e) e)
    right <- tryCatch(do.call(fit_changed, spec$right), error = function(e) e)
    if (inherits(left, "error") || inherits(right, "error")) {
      failure <- if (inherits(left, "error")) conditionMessage(left) else conditionMessage(right)
      return(tibble::tibble(
        check = spec$check,
        comparison = "nested_model",
        passed = FALSE,
        max_abs_beta_change = NA_real_,
        max_abs_IIE_change = NA_real_,
        max_abs_total_IIE_change = NA_real_,
        max_abs_comparison_difference = NA_real_,
        failure_reason = failure
      ))
    }

    max_diff <- sensitivity_B_compare_changed(left, right)
    tibble::tibble(
      check = spec$check,
      comparison = "nested_model",
      passed = max_diff <= tolerance,
      max_abs_beta_change = NA_real_,
      max_abs_IIE_change = NA_real_,
      max_abs_total_IIE_change = NA_real_,
      max_abs_comparison_difference = max_diff,
      failure_reason = NA_character_
    )
  })

  dplyr::bind_rows(baseline_result, pair_result)
}

sensitivity_B_representative_B3 <- function(data,
                                            baseline,
                                            p_U = 0.3,
                                            scenarios = tibble::tribble(
                                              ~s_M, ~s_Y,
                                              -0.5, -0.5,
                                              0.5, 0.5,
                                              0.5, -0.5
                                            ),
                                            maxit = 1000,
                                            reltol = 1e-8) {
  purrr::pmap_dfr(scenarios, function(s_M, s_Y) {
    fit <- fit_sensitivity_B_one(
      data = data,
      baseline = baseline,
      scenario = "B3_MY",
      p_U = p_U,
      OR_U_adj = 1,
      OR_U_far = 1,
      s_M = s_M,
      s_Y = s_Y,
      maxit = maxit,
      reltol = reltol
    )
    changed <- sensitivity_B_add_changes(fit$summary, fit$path, baseline)
    changed$path |>
      dplyr::mutate(
        max_abs_beta_change_scenario = max(abs(beta_change), na.rm = TRUE),
        max_abs_IIE_change_scenario = max(abs(IIE_change), na.rm = TRUE)
      )
  })
}

build_sensitivity_B_plots <- function(result_list) {
  make_heatmap <- function(df, title) {
    if (nrow(df) == 0 || !"s_M1" %in% names(df)) {
      return(NULL)
    }
    ggplot2::ggplot(
      df,
      ggplot2::aes(x = s_M1, y = s_Y, fill = total_IIE_change)
    ) +
      ggplot2::geom_tile() +
      ggplot2::facet_grid(contrast ~ OR_U_adj) +
      ggplot2::scale_fill_gradient2() +
      ggplot2::labs(x = "s_M", y = "s_Y", fill = "Change in total IIE", title = title) +
      ggplot2::theme_bw(base_size = 12)
  }

  purrr::imap(result_list, function(res, name) {
    make_heatmap(res$summary, name)
  })
}

run_sensitivity_B_latent_U <- function(data,
                                       p_U_grid = 0.3,
                                       OR_UX_grid = c(1, 1.5, 2, 3),
                                       s_M_grid = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5),
                                       s_Y_grid = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5),
                                       scenarios = c("B1_XY", "B2_XM", "B3_MY", "B4_global"),
                                       maxit = 1000,
                                       reltol = 1e-8,
                                       null_tolerance = 1e-5) {
  check_sensitivity_B_packages()
  baseline <- sensitivity_B_baseline_fit(data)

  total_time <- system.time({
    null_validation <- validate_sensitivity_B_null(
      data = data,
      baseline = baseline,
      tolerance = null_tolerance,
      maxit = maxit,
      reltol = reltol
    )
    scenario_null_validation <- validate_sensitivity_B_scenario_nulls(
      data = data,
      baseline = baseline,
      tolerance = max(null_tolerance * 10, 1e-4),
      maxit = maxit,
      reltol = reltol
    )
    nesting_validation <- validate_sensitivity_B_nesting(
      data = data,
      baseline = baseline,
      tolerance = max(null_tolerance * 10, 1e-4),
      nonzero_tolerance = max(null_tolerance * 10, 1e-4),
      maxit = maxit,
      reltol = reltol
    )
    representative_B3 <- sensitivity_B_representative_B3(
      data = data,
      baseline = baseline,
      maxit = maxit,
      reltol = reltol
    )

    scenario_results <- purrr::map(scenarios, function(scenario) {
      grid <- make_sensitivity_B_grid(
        scenario = scenario,
        p_U_grid = p_U_grid,
        OR_UX_grid = OR_UX_grid,
        s_M_grid = s_M_grid,
        s_Y_grid = s_Y_grid
      )
      run_sensitivity_B_grid(
        data = data,
        baseline = baseline,
        scenario = scenario,
        grid = grid,
        maxit = maxit,
        reltol = reltol
      )
    }) |>
      stats::setNames(scenarios)

    decomposition_validation <- validate_sensitivity_B_decomposition_identity(
      scenario_results,
      tolerance = max(null_tolerance, 1e-8)
    )
    path_behavior_validation <- validate_sensitivity_B_path_behavior(
      scenario_results,
      tolerance = max(null_tolerance * 10, 1e-4),
      nonzero_tolerance = max(null_tolerance * 10, 1e-4)
    )
  })

  timing <- purrr::imap_dfr(scenario_results, function(res, scenario) {
    tibble::tibble(
      scenario = scenario,
      grid_points = nrow(res$grid),
      joint_elapsed_seconds = sum(res$convergence$joint_elapsed_seconds, na.rm = TRUE),
      total_elapsed_seconds = sum(res$convergence$total_elapsed_seconds, na.rm = TRUE)
    )
  }) |>
    dplyr::bind_rows(tibble::tibble(
      scenario = "all",
      grid_points = sum(vapply(scenario_results, function(x) nrow(x$grid), integer(1))),
      joint_elapsed_seconds = NA_real_,
      total_elapsed_seconds = unname(total_time[["elapsed"]])
    ))

  list(
    baseline = baseline,
    sensitivity_B1_XY = scenario_results$B1_XY,
    sensitivity_B2_XM = scenario_results$B2_XM,
    sensitivity_B3_MY = scenario_results$B3_MY,
    sensitivity_B4_global = scenario_results$B4_global,
    null_validation = null_validation,
    scenario_null_validation = scenario_null_validation,
    nesting_validation = nesting_validation,
    decomposition_validation = decomposition_validation,
    path_behavior_validation = path_behavior_validation,
    representative_B3 = representative_B3,
    plots = build_sensitivity_B_plots(scenario_results),
    timing = timing,
    likelihood_specification = paste(
      "For fixed p_U, OR_U_adj, OR_U_far, lambda_U, and delta_U,",
      "the estimator maximizes sum_i log[sum_{u=0}^1 P(U_i=u|X_i)",
      "f_3(M_i|X_i,C_i,U_i=u; Sigma_M) f(Y_i|X_i,M_i,C_i,U_i=u; sigma_Y^2)],",
      "with component means M = X_M alpha - lambda_U U and Y = X_Y theta - delta_U U.",
      "The primary reported mediation TE is TE = direct_effect + total_IIE;",
      "the separately fitted reduced-form total-model coefficient is retained as TE_reduced_form.",
      "Optimization uses stats::optim(method='BFGS') on regression coefficients,",
      "Cholesky-parameterized Sigma_M, and log residual SDs."
    )
  )
}
