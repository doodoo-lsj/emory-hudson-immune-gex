suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

find_project_root <- function() {
  candidates <- unique(normalizePath(
    c(getwd(), file.path(getwd(), "..")),
    winslash = "/",
    mustWork = FALSE
  ))

	  for (candidate in candidates) {
	    if (file.exists(file.path(candidate, "R", "pca_mediation_pipeline.R")) &&
	        file.exists(file.path(candidate, "R", "data_input_app.R")) &&
	        file.exists(file.path(candidate, "R", "bayesian_mediation_v1.R")) &&
	        file.exists(file.path(candidate, "R", "bayesian_sensitivity_A_only.R")) &&
	        file.exists(file.path(candidate, "R", "analysis_bayesian_app.R")) &&
	        file.exists(file.path(candidate, "R", "bayesian_dynamic_app.R")) &&
	        file.exists(file.path(candidate, "R", "frequentist_dynamic_app.R")) &&
	        file.exists(file.path(candidate, "R", "sensitivity_v2.R")) &&
	        file.exists(file.path(candidate, "data", "raw", "pt16_emory_GEX_immune_FULL_v2.rds")) &&
	        file.exists(file.path(candidate, "results", "bayesian", "current", "emory_bayesian_v1_sensitivity_A_app_artifact.rds"))) {
	      return(candidate)
	    }
	  }

  stop("Could not locate the project root from the Shiny app directory.")
}

project_root <- find_project_root()
source(file.path(project_root, "R", "pca_mediation_pipeline.R"))
source(file.path(project_root, "R", "data_input_app.R"))
source(file.path(project_root, "R", "bayesian_mediation_v1.R"))
source(file.path(project_root, "R", "bayesian_sensitivity_A_only.R"))
source(file.path(project_root, "R", "analysis_bayesian_app.R"))
source(file.path(project_root, "R", "bayesian_dynamic_app.R"))
source(file.path(project_root, "R", "frequentist_dynamic_app.R"))
source(file.path(project_root, "R", "sensitivity_v2.R"))

rds_path <- file.path(project_root, "data", "raw", "pt16_emory_GEX_immune_FULL_v2.rds")
bayesian_artifact_path <- default_bayesian_app_artifact_path(project_root)
bayesian_app_artifact <- load_bayesian_app_artifact(bayesian_artifact_path)
example_inputs <- example_data_inputs(rds_path)
example_spec <- build_analysis_spec(
  data_source = "example",
  main_data = example_inputs$main_data,
  feature_data = example_inputs$feature_data,
  roles = example_inputs$defaults,
  feature_id = example_inputs$defaults$feature_observation_id,
  excluded_feature_columns = example_inputs$defaults$excluded_feature_columns
)
example_analysis <- compute_analysis_from_spec(example_spec)
coord_x <- "imagecol"
coord_y <- "imagerow"
exposures <- c("X_adjacent", "X_far")
mediators <- c("PC1_R", "PC2_R", "PC3")
covariates <- c("x_coord", "y_coord")
rho_grid <- seq(-0.5, 0.5, by = 0.05)
bayesian_sensitivity_rho_grid <- seq(-0.3, 0.3, by = 0.05)
validated_pc_ids <- c("PC1", "PC2", "PC3")
pc_display_names <- c(PC1_R = "PC1", PC2_R = "PC2", PC3 = "PC3")
pc_mediator_ids <- c(PC1 = "PC1_R", PC2 = "PC2_R", PC3 = "PC3")

contrast_label <- function(x) {
  dplyr::case_when(
    x == "X_adjacent" ~ "adjacent vs inside",
    x == "X_far" ~ "far vs inside",
    TRUE ~ x
  )
}

round_numeric <- function(df, digits = 4) {
  df |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, digits)))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

mediator_label <- function(x) {
  out <- unname(pc_display_names[x])
  out[x == "all"] <- "All mediators"
  out[grepl("^PC[0-9]+$", x)] <- x[grepl("^PC[0-9]+$", x)]
  ifelse(is.na(out), x, out)
}

mediator_filter_choices <- function(mediator_ids = mediators) {
  c("All mediators" = "all", stats::setNames(mediator_ids, mediator_label(mediator_ids)))
}

contrast_filter_choices <- function() {
  c("All contrasts" = "all", "adjacent vs inside" = "adjacent vs inside", "far vs inside" = "far vs inside")
}

pc_sign_vector <- function(input) {
  stats::setNames(
    vapply(seq_len(6), function(i) {
      value <- input[[paste0("pc", i, "_orientation")]]
      if (identical(value, "reverse")) -1 else 1
    }, numeric(1)),
    paste0("PC", seq_len(6))
  )
}

mediation_sign_vector <- function(input) {
  pc_signs <- pc_sign_vector(input)
  c(
    PC1_R = pc_signs[["PC1"]],
    PC2_R = pc_signs[["PC2"]],
    PC3 = pc_signs[["PC3"]]
  )
}

pc_candidate_ids <- function(analysis_obj, max_pcs = 6) {
  paste0("PC", seq_len(min(max_pcs, ncol(analysis_obj$pca$pca_fit$x))))
}

selected_pc_ids <- function(input, analysis_obj = NULL) {
  candidates <- if (is.null(analysis_obj)) {
    paste0("PC", seq_len(6))
  } else {
    pc_candidate_ids(analysis_obj)
  }

  selected <- NULL
  pc_flags <- vapply(seq_along(candidates), function(i) {
    isTRUE(input[[paste0("pca_pc", i, "_selected")]])
  }, logical(1))
  if (any(!vapply(seq_along(candidates), function(i) {
    is.null(input[[paste0("pca_pc", i, "_selected")]])
  }, logical(1)))) {
    selected <- candidates[pc_flags]
  }
  selected <- selected %||% input$pca_mediation_pcs %||% candidates
  if (!is.null(analysis_obj)) {
    selected <- intersect(selected, candidates)
  }
  selected
}

validated_selection_selected <- function(selected) {
  setequal(selected, validated_pc_ids)
}

orient_spatial_scores <- function(plot_df, signs, pca_fit = NULL) {
  out <- plot_df |>
    mutate(
      PC1 = PC1_R * signs[["PC1"]],
      PC2 = PC2_R * signs[["PC2"]],
      PC3 = PC3 * signs[["PC3"]]
    )

  if (!is.null(pca_fit)) {
    for (i in 4:min(6, ncol(pca_fit$x))) {
      pc <- paste0("PC", i)
      out[[pc]] <- pca_fit$x[, i] * signs[[pc]]
    }
  }

  out
}

orient_loading_data <- function(loading_df, signs, pca_fit = NULL) {
  out <- loading_df |>
    mutate(
      PC1 = PC1_R * signs[["PC1"]],
      PC2 = PC2_R * signs[["PC2"]],
      PC3 = PC3 * signs[["PC3"]]
    )

  if (!is.null(pca_fit)) {
    rotation <- pca_fit$rotation
    for (i in 4:min(6, ncol(rotation))) {
      pc <- paste0("PC", i)
      out[[pc]] <- rotation[, i] * signs[[pc]]
    }
  }

  out
}

orient_interval_columns <- function(df, prefix, signs) {
  sign_vec <- signs[df$mediator]
  old_lwr <- df[[paste0(prefix, "_lwr")]]
  old_upr <- df[[paste0(prefix, "_upr")]]
  reversed <- sign_vec < 0

  df[[paste0(prefix, "_est")]] <- df[[paste0(prefix, "_est")]] * sign_vec
  df[[paste0(prefix, "_lwr")]] <- ifelse(reversed, -old_upr, old_lwr)
  df[[paste0(prefix, "_upr")]] <- ifelse(reversed, -old_lwr, old_upr)
  df
}

orient_bayesian_path <- function(df, signs) {
  if (!all(df$mediator %in% names(signs))) {
    df$mediator_id <- df$mediator
    df$mediator <- mediator_label(df$mediator)
    return(df)
  }
  df <- orient_interval_columns(df, "alpha", signs)
  df <- orient_interval_columns(df, "beta", signs)
  df$mediator_id <- df$mediator
  df$mediator <- mediator_label(df$mediator)
  df
}

orient_frequentist_path <- function(df, signs) {
  sign_vec <- signs[df$mediator]
  df |>
    mutate(
      mediator_id = mediator,
      mediator = mediator_label(mediator),
      alpha = alpha * sign_vec,
      beta = beta * sign_vec
    )
}

filter_by_mediation_inputs <- function(df, contrast = "all", mediator = "all") {
  out <- df
  if (!identical(contrast, "all") && "contrast" %in% names(out)) {
    out <- dplyr::filter(out, contrast == !!contrast)
  }
  if (!identical(mediator, "all") && "mediator_id" %in% names(out)) {
    out <- dplyr::filter(out, mediator_id == !!mediator)
  }
  out
}

format_number <- function(x, digits = 3) {
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x)
  out[ok] <- formatC(x[ok], format = "f", digits = digits)
  out
}

format_estimate_interval <- function(est, lwr = NULL, upr = NULL, digits = 3) {
  est_txt <- format_number(est, digits)
  if (is.null(lwr) || is.null(upr)) {
    return(est_txt)
  }

  paste0(est_txt, " [", format_number(lwr, digits), ", ", format_number(upr, digits), "]")
}

present_bayesian_mediation_summary <- function(df) {
  df |>
    transmute(
      Contrast = contrast,
      `Total Effect` = format_estimate_interval(TE_est, TE_lwr, TE_upr),
      `Direct Effect` = format_estimate_interval(direct_est, direct_lwr, direct_upr),
      `Indirect Effect (Total)` = format_estimate_interval(total_IIE_est, total_IIE_lwr, total_IIE_upr),
      `Proportion Mediated` = format_estimate_interval(PM_est, PM_lwr, PM_upr)
    )
}

present_bayesian_mediation_probabilities <- function(df) {
  out <- tibble::tibble(Contrast = df$contrast)
  if ("TE_Pr_gt_0" %in% names(df)) {
    out$`Total Effect` <- format_number(df$TE_Pr_gt_0)
  }
  if ("direct_Pr_gt_0" %in% names(df)) {
    out$`Direct Effect` <- format_number(df$direct_Pr_gt_0)
  }
  if ("total_IIE_Pr_gt_0" %in% names(df)) {
    out$`Indirect Effect` <- format_number(df$total_IIE_Pr_gt_0)
  }

  out
}

present_bayesian_mediation_path <- function(df) {
  out <- df |>
    transmute(
      Contrast = contrast,
      Mediator = mediator,
      `Exposure -> Mediator` = format_estimate_interval(alpha_est, alpha_lwr, alpha_upr),
      `Mediator -> Outcome` = format_estimate_interval(beta_est, beta_lwr, beta_upr),
      `Indirect Effect` = format_estimate_interval(IIE_est, IIE_lwr, IIE_upr)
    )

  if ("IIE_Pr_gt_0" %in% names(df)) {
    out$`Pr(Indirect Effect > 0)` <- format_number(df$IIE_Pr_gt_0)
  }
  if ("IIE_Pr_lt_0" %in% names(df)) {
    out$`Pr(Indirect Effect < 0)` <- format_number(df$IIE_Pr_lt_0)
  }

  out
}

present_frequentist_mediation_summary <- function(df) {
  df |>
    transmute(
      Contrast = contrast,
      `Total Effect` = format_estimate_interval(TE),
      `Direct Effect` = format_estimate_interval(direct_effect),
      `Indirect Effect (Total)` = format_estimate_interval(total_indirect_effect),
      `Proportion Mediated` = format_estimate_interval(PM)
    )
}

present_frequentist_mediation_path <- function(df) {
  df |>
    transmute(
      Contrast = contrast,
      Mediator = mediator,
      `Exposure -> Mediator` = format_estimate_interval(alpha, digits = 4),
      `Mediator -> Outcome` = format_estimate_interval(beta, digits = 4),
      `Indirect Effect` = format_estimate_interval(indirect_component, digits = 4)
    )
}

present_bootstrap_summary_ci <- function(df, contrast_labels = NULL) {
  label <- if ("contrast" %in% names(df)) {
    df$contrast
  } else if (is.null(contrast_labels)) {
    contrast_label(df$exposure)
  } else {
    out <- unname(contrast_labels[df$exposure])
    ifelse(is.na(out), df$exposure, out)
  }

  df |>
    mutate(contrast = label) |>
    transmute(
      Contrast = contrast,
      `Total Effect` = format_estimate_interval(TE_mean, TE_lwr, TE_upr),
      `Direct Effect` = format_estimate_interval(NDE_mean, NDE_lwr, NDE_upr),
      `Indirect Effect (Total)` = format_estimate_interval(NIE_mean, NIE_lwr, NIE_upr),
      `Proportion Mediated` = format_estimate_interval(PM_mean, PM_lwr, PM_upr)
    )
}

present_bootstrap_path_ci <- function(df, contrast_labels = NULL) {
  label <- if ("contrast" %in% names(df)) {
    df$contrast
  } else if (is.null(contrast_labels)) {
    contrast_label(df$exposure)
  } else {
    out <- unname(contrast_labels[df$exposure])
    ifelse(is.na(out), df$exposure, out)
  }

  df |>
    mutate(contrast = label) |>
    transmute(
      Contrast = contrast,
      Mediator = mediator_label(mediator),
      `Exposure -> Mediator` = format_estimate_interval(alpha_mean, alpha_lwr, alpha_upr),
      `Mediator -> Outcome` = format_estimate_interval(beta_mean, beta_lwr, beta_upr),
      `Indirect Effect` = format_estimate_interval(indirect_mean, indirect_lwr, indirect_upr)
    )
}

present_bayesian_sensitivity_summary <- function(df) {
  df |>
    transmute(
      Contrast = contrast,
      `Direct Effect` = format_estimate_interval(direct_est, direct_lwr, direct_upr),
      `Indirect Effect (Total)` = format_estimate_interval(total_IIE_est, total_IIE_lwr, total_IIE_upr),
      `Total Effect` = format_estimate_interval(TE_est, TE_lwr, TE_upr),
      `Proportion Mediated` = format_estimate_interval(PM_est, PM_lwr, PM_upr)
    )
}

present_bayesian_sensitivity_path <- function(df) {
  df |>
    transmute(
      Contrast = contrast,
      rho = format_number(rho, digits = 2),
      `Mediator Varied` = mediator_label(mediator_varied),
      Mediator = mediator_label(mediator),
      `Indirect Effect` = format_estimate_interval(IIE_est, IIE_lwr, IIE_upr),
      `Pr(Indirect Effect > 0)` = format_number(IIE_Pr_gt_0),
      `Pr(Indirect Effect < 0)` = format_number(IIE_Pr_lt_0)
    )
}

present_frequentist_sensitivity_selected <- function(df) {
  if ("varied_mediator" %in% names(df)) {
    df |>
      transmute(
        Contrast = contrast,
        `Direct Effect` = format_estimate_interval(adjusted_direct_effect),
        `Indirect Effect (Total)` = format_estimate_interval(adjusted_total_indirect),
        `Total Effect` = format_estimate_interval(adjusted_total_effect),
        `Proportion Mediated` = format_estimate_interval(PM_adjusted)
      )
  } else {
    df |>
      transmute(
        Contrast = contrast,
        `Direct Effect` = format_estimate_interval(adjusted_direct_effect),
        `Indirect Effect (Total)` = format_estimate_interval(adjusted_total_indirect),
        `Total Effect` = format_estimate_interval(adjusted_total_effect),
        `Proportion Mediated` = format_estimate_interval(PM_adjusted)
      )
  }
}

present_feasibility_status <- function(summary_df) {
  fraction <- min(summary_df$valid_draw_fraction, na.rm = TRUE)
  percent <- round(100 * fraction)

  if (!is.finite(fraction) || fraction <= 0) {
    list(status = "Infeasible", text = "Infeasible for the selected rho. No posterior draws support this setting.")
  } else if (fraction < 1) {
    list(
      status = "Partially feasible",
      text = paste0("Feasible for ", percent, "% of posterior draws. Interpret this rho setting with caution.")
    )
  } else {
    list(status = "Feasible", text = "Feasible for 100% of posterior draws.")
  }
}

pm_warning_needed <- function(summary_df) {
  if (nrow(summary_df) == 0) {
    return(FALSE)
  }

  te_near_zero <- abs(summary_df$TE_est) < 0.01
  te_crosses_zero <- summary_df$TE_lwr <= 0 & summary_df$TE_upr >= 0
  opposing_signs <- summary_df$direct_est * summary_df$total_IIE_est < 0

  any(te_near_zero | te_crosses_zero | opposing_signs, na.rm = TRUE)
}

pm_warning_needed_frequentist <- function(summary_df) {
  if (nrow(summary_df) == 0) {
    return(FALSE)
  }
  te_near_zero <- abs(summary_df$adjusted_total_effect) < 0.01
  opposing_signs <- summary_df$adjusted_direct_effect * summary_df$adjusted_total_indirect < 0
  any(te_near_zero | opposing_signs, na.rm = TRUE)
}

spatial_plot_limits <- function(data, padding_fraction = 0.06) {
  x_range <- range(data$x_plot, na.rm = TRUE)
  y_range <- range(data$y_plot, na.rm = TRUE)
  x_pad <- diff(x_range) * padding_fraction
  y_pad <- diff(y_range) * padding_fraction
  if (!is.finite(x_pad) || x_pad == 0) x_pad <- 1
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1

  list(
    x = c(x_range[[1]] - x_pad, x_range[[2]] + x_pad),
    y = c(y_range[[1]] - y_pad, y_range[[2]] + y_pad)
  )
}

compute_current_analysis <- function() {
  example_analysis
}

compute_sensitivity <- function(dat_full_allpc,
                                exposures = exposures,
                                mediators = mediators,
                                covariates = covariates) {
  common <- residual_corr_sensitivity(
    data = dat_full_allpc,
    exposures = exposures,
    outcome = "Y",
    mediators = mediators,
    covariates = covariates,
    rho_grid = rho_grid,
    mode = "common"
  )

  one_at_a_time <- residual_corr_sensitivity(
    data = dat_full_allpc,
    exposures = exposures,
    outcome = "Y",
    mediators = mediators,
    covariates = covariates,
    rho_grid = rho_grid,
    mode = "one_at_a_time"
  )

  list(
    common = common,
    one_at_a_time = one_at_a_time,
    rho0_common = validate_rho0_matches_observed(common),
    rho0_one_at_a_time = validate_rho0_matches_observed(one_at_a_time),
    tipping_common = summarize_rho_tipping_points(common),
    tipping_one_at_a_time = summarize_rho_tipping_points(one_at_a_time),
    plots = build_rho_sensitivity_plots(common, one_at_a_time)
  )
}

analysis_cache <- local({
  value <- example_analysis
  function() {
    value
  }
})

sensitivity_cache <- local({
  value <- NULL
  function(dat_full_allpc) {
    if (is.null(value)) {
      value <<- compute_sensitivity(dat_full_allpc)
    }
    value
  }
})

compute_bayesian_sensitivity_curve <- function(artifact, mode = c("common", "one_at_a_time"), mediator = NULL) {
  mode <- match.arg(mode)
  valid_mediators <- bayes_A_bundle_mediators(artifact$sensitivity_A$residual_draws)
  raw <- compute_bayesian_sensitivity_A_from_draws(
    residual_draws = artifact$sensitivity_A$residual_draws,
    rho_grid = bayesian_sensitivity_rho_grid,
    mode = mode,
    mediators = valid_mediators,
    mediators_varied = if (identical(mode, "one_at_a_time") && !is.null(mediator)) mediator else valid_mediators
  )

  summary <- canonicalize_bayesian_A_summary(raw)
  path <- canonicalize_bayesian_A_path(raw)
  feasibility <- raw$feasibility

  if (identical(mode, "one_at_a_time")) {
    mediator <- mediator %||% valid_mediators[[1]]
    summary <- dplyr::filter(summary, mediator_varied == mediator)
    path <- dplyr::filter(path, mediator_varied == mediator)
    feasibility <- dplyr::filter(feasibility, mediator_varied == mediator)
  }

  list(summary = summary, path = path, feasibility = feasibility)
}

bayesian_sensitivity_curve_cache <- local({
  values <- list()
  function(artifact, artifact_key, mode = c("common", "one_at_a_time"), mediator = NULL) {
    mode <- match.arg(mode)
    key <- paste(artifact_key, mode, mediator %||% "all", sep = ":")
    if (is.null(values[[key]])) {
      values[[key]] <<- compute_bayesian_sensitivity_curve(artifact = artifact, mode = mode, mediator = mediator)
    }
    values[[key]]
  }
})

bayesian_sensitivity_all_mediator_curve_cache <- local({
  values <- list()
  function(artifact, artifact_key) {
    key <- paste(artifact_key, "one_at_a_time", "all_mediators", sep = ":")
    if (is.null(values[[key]])) {
      raw <- compute_bayesian_sensitivity_A_from_draws(
        residual_draws = artifact$sensitivity_A$residual_draws,
        rho_grid = bayesian_sensitivity_rho_grid,
        mode = "one_at_a_time",
        mediators = bayes_A_bundle_mediators(artifact$sensitivity_A$residual_draws),
        mediators_varied = bayes_A_bundle_mediators(artifact$sensitivity_A$residual_draws)
      )
      values[[key]] <<- list(
        summary = canonicalize_bayesian_A_summary(raw),
        path = canonicalize_bayesian_A_path(raw),
        feasibility = raw$feasibility
      )
    }
    values[[key]]
  }
})

bayesian_sensitivity_selected_cache <- local({
  values <- list()
  function(artifact, artifact_key, rho, mode = c("common", "one_at_a_time"), mediator = NULL) {
    mode <- match.arg(mode)
    key <- paste(artifact_key, mode, mediator %||% "all", sprintf("%.12g", rho), sep = ":")
    if (is.null(values[[key]])) {
      values[[key]] <<- compute_bayesian_app_sensitivity_A(
        artifact = artifact,
        rho = rho,
        mode = mode,
        mediator = mediator
      )
    }
    values[[key]]
  }
})

top_loading_table <- function(loading_df, pc, direction = c("positive", "negative"), n = 30) {
  direction <- match.arg(direction)
  if (direction == "positive") {
    loading_df |>
      arrange(desc(.data[[pc]])) |>
      slice_head(n = n) |>
      select(gene, loading = all_of(pc))
  } else {
    loading_df |>
      arrange(.data[[pc]]) |>
      slice_head(n = n) |>
      select(gene, loading = all_of(pc))
  }
}

app_panel_card <- function(title, ..., class = "") {
  tags$section(
    class = paste("app-card", class),
    h4(title),
    ...
  )
}

app_tutorial_section <- function(title, ..., open = FALSE) {
  args <- c(
    list(class = "tutorial-section"),
    if (isTRUE(open)) list(open = NA) else list(),
    list(tags$summary(title), tags$div(class = "tutorial-body", ...))
  )
  do.call(tags$details, args)
}

app_tutorial_step <- function(number, title, text) {
  tags$div(
    class = "tutorial-step",
    tags$div(class = "tutorial-step-number", number),
    tags$div(
      class = "tutorial-step-copy",
      tags$h4(title),
      tags$p(text)
    )
  )
}

app_tutorial_tile <- function(title, text) {
  tags$div(
    class = "tutorial-tile",
    tags$h4(title),
    tags$p(text)
  )
}

app_setup_section <- function(title, ...) {
  tags$section(
    class = "setup-sidebar-section",
    h4(title),
    ...
  )
}

app_kv_item <- function(label, value, class = "") {
  tags$div(
    class = paste("app-kv-item", class),
    tags$span(class = "app-kv-label", label),
    tags$span(class = "app-kv-value", value)
  )
}

pc_orientation_label <- function(sign_value) {
  if (isTRUE(sign_value < 0)) "Reverse" else "Keep"
}

pca_spec_label <- function(spec) {
  if (is.null(spec) || length(spec$selected_pcs) == 0) {
    return("Not confirmed")
  }
  paste(spec$selected_pcs, collapse = ", ")
}

pca_orientation_summary <- function(selected_pcs, pc_signs) {
  if (length(selected_pcs) == 0) {
    return("None")
  }
  paste(
    vapply(
      selected_pcs,
      function(pc) paste0(pc, ": ", pc_orientation_label(pc_signs[[pc]])),
      character(1)
    ),
    collapse = "; "
  )
}

ui <- navbarPage(
  title = "PCA Mediation v0",
  id = "main_nav",
  header = tagList(
    tags$head(
      tags$script(HTML(
        "
        Shiny.addCustomMessageHandler('setDisabled', function(x) {
          $('#' + x.id).prop('disabled', x.disabled);
        });
        Shiny.addCustomMessageHandler('setNavWorkflow', function(x) {
          var navItems = $('a[data-value]');
          navItems
            .removeClass('nav-complete nav-available nav-unavailable nav-utility')
            .removeAttr('aria-disabled')
            .removeAttr('tabindex');
          navItems.find('.nav-check').remove();

          x.steps.forEach(function(step) {
            var anchor = $('a[data-value=\"' + step.value + '\"]').first();
            if (anchor.length === 0) return;
            anchor.addClass(step.workflow ? 'nav-available' : 'nav-utility');
            if (step.completed) {
              anchor.addClass('nav-complete');
              anchor.prepend('<span class=\"nav-check\">&#10003;</span>');
            }
            if (!step.available) {
              anchor.addClass('nav-unavailable');
              anchor.attr('aria-disabled', 'true');
              anchor.attr('tabindex', '-1');
            }
          });
        });
        "
      )),
      tags$style(HTML("
	      :root {
	        --app-ink: #243447;
		        --app-accent: #274690;
		        --app-accent-2: #2f80ed;
		        --app-accent-soft: #eef4ff;
		        --app-blue: #2f80ed;
		        --app-blue-soft: #eaf2ff;
		        --app-violet: #7c5cc4;
		        --app-violet-soft: #f1edff;
		        --app-red: #c2413a;
		        --app-red-soft: #fdeceb;
		        --app-slate-accent: #64748b;
		        --app-slate-soft: #f1f5f9;
		        --app-border: #dde8ef;
		        --app-muted: #5f6f7f;
		        --app-soft: #f7fafc;
		      }
	      body {
	        background: #f7fafc;
	        color: var(--app-ink);
	      }
	      .navbar {
	        margin-bottom: 0;
	        border-bottom: 1px solid var(--app-border);
		        box-shadow: 0 2px 8px rgba(39, 70, 144, 0.08);
	      }
	      .navbar-default {
	        background: #ffffff;
	        border-color: transparent transparent var(--app-border) transparent;
	      }
	      .navbar-default .navbar-brand {
	        color: var(--app-accent);
	        font-weight: 800;
	      }
	      .navbar-default .navbar-brand:hover,
	      .navbar-default .navbar-brand:focus {
		        color: #1d3470;
	      }
	      .navbar-default .navbar-nav > li > a {
	        color: var(--app-ink);
	        font-weight: 600;
	        border-bottom: 3px solid transparent;
	        transition: background-color 120ms ease, color 120ms ease, border-color 120ms ease;
	      }
	      .navbar-default .navbar-nav > li.active > a,
	      .navbar-default .navbar-nav > li.active > a:focus,
	      .navbar-default .navbar-nav > li.active > a:hover {
	        background: var(--app-accent-soft);
	        color: var(--app-accent);
	        border-bottom-color: var(--app-accent-2);
	      }
      .navbar-default .navbar-nav > li > a.nav-unavailable {
        color: #9aa7b3;
        cursor: default;
        pointer-events: none;
      }
      .navbar-default .navbar-nav > li > a.nav-unavailable:hover,
      .navbar-default .navbar-nav > li > a.nav-unavailable:focus {
        background: transparent;
        color: #9aa7b3;
      }
	      .nav-check {
		        color: var(--app-blue);
	        font-size: 11px;
	        margin-right: 5px;
        position: relative;
        top: -1px;
      }
      .setup-shell {
        display: grid;
        grid-template-columns: minmax(280px, 340px) minmax(0, 1fr);
        gap: 18px;
        align-items: start;
        margin-top: 12px;
      }
	      .setup-sidebar {
	        position: sticky;
	        top: 12px;
	        background: #ffffff;
	        border: 1px solid var(--app-border);
	        border-top: 3px solid var(--app-accent);
	        border-radius: 8px;
	        padding: 14px;
		        box-shadow: 0 4px 12px rgba(39, 70, 144, 0.07);
	      }
      .setup-main {
        min-width: 0;
      }
      .setup-sidebar-section {
        border-bottom: 1px solid #e8eef5;
        padding-bottom: 12px;
        margin-bottom: 14px;
      }
      .setup-sidebar-section:last-child {
        border-bottom: 0;
        margin-bottom: 0;
        padding-bottom: 0;
      }
	      .setup-sidebar-section h4,
	      .app-card h4 {
	        margin-top: 0;
	        margin-bottom: 10px;
	        font-size: 15px;
	        font-weight: 700;
	        color: var(--app-ink);
	      }
      .setup-sidebar .form-group {
        margin-bottom: 10px;
      }
	      .setup-sidebar .btn-primary {
	        width: 100%;
	      }
		      .btn-primary,
		      .btn.btn-primary,
		      .action-button.btn-primary,
		      input.btn-primary {
		        background-color: var(--app-accent) !important;
		        border-color: var(--app-accent) !important;
		        color: #ffffff !important;
		        font-weight: 700;
		      }
		      .btn-primary:visited,
		      .btn-primary:active,
		      .btn.btn-primary:visited,
		      .btn.btn-primary:active,
		      .action-button.btn-primary:visited,
		      .action-button.btn-primary:active {
		        color: #ffffff !important;
		      }
		      .btn-primary:hover,
		      .btn-primary:focus,
		      .btn.btn-primary:hover,
		      .btn.btn-primary:focus,
		      .action-button.btn-primary:hover,
		      .action-button.btn-primary:focus,
		      input.btn-primary:hover,
		      input.btn-primary:focus {
		        background-color: #1f3978 !important;
		        border-color: #1f3978 !important;
		        color: #ffffff !important;
		      }
		      .btn-primary[disabled],
		      .btn-primary.disabled,
		      .btn.btn-primary[disabled],
		      .btn.btn-primary.disabled,
		      .action-button.btn-primary[disabled],
		      .action-button.btn-primary.disabled,
		      input.btn-primary[disabled],
		      input.btn-primary.disabled {
		        background-color: #6f7fb3 !important;
		        border-color: #6f7fb3 !important;
		        color: #ffffff !important;
		        opacity: 0.88;
		      }
	      .btn-default {
	        background: #ffffff;
		        border-color: #c6d6f5;
	        color: var(--app-accent);
	        font-weight: 700;
	      }
	      .btn-default:hover,
	      .btn-default:focus {
	        background: var(--app-accent-soft);
		        border-color: var(--app-blue);
	        color: var(--app-accent);
	      }
		      .app-card {
		        background: #ffffff;
		        border: 1px solid var(--app-border);
		        border-radius: 8px;
	        padding: 14px;
	        margin-bottom: 14px;
		        box-shadow: 0 2px 8px rgba(22, 34, 51, 0.04);
		        min-width: 0;
		      }
		      .tutorial-shell {
		        max-width: 1080px;
		        margin: 18px auto 0 auto;
		      }
			      .tutorial-hero {
				        background: #f6f8ff;
				        border: 1px solid #d9e3ff;
				        border-left: 5px solid var(--app-blue);
			        border-radius: 8px;
			        padding: 22px 24px;
			        margin-bottom: 14px;
				        box-shadow: 0 5px 18px rgba(39, 70, 144, 0.10);
			      }
			      .tutorial-eyebrow {
				        color: var(--app-accent);
			        font-size: 12px;
		        font-weight: 700;
		        text-transform: uppercase;
		        margin-bottom: 6px;
		      }
			      .tutorial-hero h2 {
			        margin: 0 0 8px 0;
			        color: var(--app-ink);
		        font-size: 28px;
		        font-weight: 700;
		        line-height: 1.15;
		      }
		      .tutorial-hero p {
		        color: var(--app-muted);
		        font-size: 15px;
		        line-height: 1.5;
		        margin: 0;
		        max-width: 780px;
		      }
	      .tutorial-workflow {
	        display: grid;
	        grid-template-columns: repeat(5, minmax(0, 1fr));
	        gap: 10px;
	        align-items: stretch;
	        margin-bottom: 14px;
	      }
	      .tutorial-step {
	        background: #ffffff;
	        border: 1px solid var(--app-border);
	        border-top: 3px solid var(--app-blue);
	        border-radius: 8px;
	        padding: 12px;
	        box-shadow: 0 2px 8px rgba(22, 34, 51, 0.04);
	        height: 100%;
	      }
			      .tutorial-step:nth-child(2),
			      .tutorial-step:nth-child(5) {
				        border-top-color: var(--app-violet);
			      }
			      .tutorial-step:nth-child(3) {
				        border-top-color: var(--app-accent);
			      }
			      .tutorial-step-number {
		        width: 26px;
		        height: 26px;
		        border-radius: 999px;
			        background: var(--app-accent-soft);
			        color: var(--app-accent);
		        display: flex;
		        align-items: center;
		        justify-content: center;
		        font-weight: 700;
		        font-size: 12px;
			        margin-bottom: 9px;
			      }
			      .tutorial-step:nth-child(2) .tutorial-step-number,
			      .tutorial-step:nth-child(5) .tutorial-step-number {
				        background: var(--app-violet-soft);
				        color: var(--app-violet);
			      }
			      .tutorial-step:nth-child(3) .tutorial-step-number {
				        background: var(--app-accent-soft);
				        color: var(--app-accent);
			      }
			      .tutorial-step h4,
			      .tutorial-tile h4 {
			        margin: 0 0 6px 0;
			        color: var(--app-ink);
		        font-size: 14px;
		        font-weight: 700;
		      }
		      .tutorial-step p,
		      .tutorial-tile p {
		        margin: 0;
		        color: var(--app-muted);
		        font-size: 12px;
		        line-height: 1.4;
		      }
	      .tutorial-tile-grid {
	        display: grid;
	        grid-template-columns: repeat(3, minmax(0, 1fr));
	        gap: 10px;
	        align-items: stretch;
	        margin-bottom: 14px;
	      }
	      .tutorial-tile {
	        background: var(--app-soft);
	        border: 1px solid #e5edf5;
	        border-top: 3px solid var(--app-blue);
	        border-radius: 8px;
	        padding: 12px;
	        height: 100%;
	      }
			      .tutorial-tile:nth-child(2) {
				        border-top-color: var(--app-accent);
				        background: var(--app-accent-soft);
			      }
			      .tutorial-tile:nth-child(3) {
				        border-top-color: var(--app-violet);
				        background: var(--app-violet-soft);
			      }
			      .tutorial-section {
			        background: #ffffff;
	        border: 1px solid var(--app-border);
	        border-radius: 8px;
	        margin-bottom: 10px;
	        box-shadow: 0 1px 2px rgba(22, 34, 51, 0.04);
	        overflow: hidden;
	      }
		      .tutorial-section summary {
		        cursor: pointer;
		        padding: 12px 14px;
		        color: var(--app-accent);
		        font-weight: 700;
		      }
	      .tutorial-body {
	        border-top: 1px solid #eef3f8;
	        padding: 12px 14px 14px 14px;
	        color: var(--app-muted);
	        font-size: 13px;
	        line-height: 1.5;
	      }
	      .tutorial-body p {
	        margin: 0 0 10px 0;
	      }
		      .tutorial-body p:last-child,
		      .tutorial-body ol,
		      .tutorial-body ul {
		        margin-bottom: 0;
		      }
	      .setup-preview-table .dataTables_wrapper {
	        width: 100%;
	      }
      .setup-preview-table .dataTables_scroll {
        overflow: visible;
      }
      .setup-preview-table .dataTables_scrollBody {
        border-bottom: 1px solid #d9e2ec;
      }
	      .setup-main-grid {
	        display: grid;
	        grid-template-columns: repeat(2, minmax(0, 1fr));
	        gap: 14px;
	        align-items: stretch;
	      }
	      .setup-main-grid > .app-card,
	      .overview-summary-grid > .app-card,
	      .pca-detail-grid > .app-card {
	        height: 100%;
	        margin-bottom: 0;
	      }
      .setup-empty {
        min-height: 220px;
        display: flex;
        flex-direction: column;
        justify-content: center;
        background: var(--app-soft);
      }
      .setup-empty p,
      .setup-hint,
      .setup-meta {
        color: var(--app-muted);
      }
      .setup-validation-list {
        list-style: none;
        padding-left: 0;
        margin-bottom: 0;
      }
      .setup-validation-list li {
        padding: 5px 0;
        border-bottom: 1px solid #eef3f8;
      }
      .setup-validation-list li:last-child {
        border-bottom: 0;
      }
      .setup-ok {
        color: #1b6e3c;
        font-weight: 700;
        margin-right: 6px;
      }
      .setup-warning {
        color: #8a5a00;
        font-weight: 700;
        margin-right: 6px;
      }
      .setup-status {
        margin-top: 8px;
        font-size: 13px;
      }
      .setup-requirement-list {
        margin: 8px 0 0 0;
        padding-left: 20px;
      }
      .setup-requirement-list > li {
        margin-bottom: 10px;
      }
      .setup-requirement-list ul {
        color: var(--app-muted);
        font-size: 13px;
        margin-top: 5px;
        padding-left: 18px;
      }
	      .workflow-context-card {
	        background: var(--app-soft);
	        border: 1px solid #e5edf5;
	        border-left: 3px solid var(--app-blue);
	        border-radius: 8px;
        padding: 11px 12px;
        margin-bottom: 10px;
      }
      .workflow-context-card h5 {
        margin: 0 0 8px 0;
        font-size: 13px;
        font-weight: 700;
        color: #243447;
      }
      .workflow-context-card .app-kv-item {
        padding: 5px 0;
        font-size: 12px;
      }
      .workflow-action {
        margin-top: 10px;
      }
      .workflow-action .btn,
      .workflow-action .btn-primary {
        width: 100%;
      }
      .workflow-note {
        color: var(--app-muted);
        font-size: 12px;
        margin: 6px 0 0 0;
      }
	      .workflow-note.warning {
	        color: #8a5a00;
	      }
	      .workflow-note.ok {
	        color: var(--app-accent);
	      }
      .overview-shell {
        width: 78%;
        max-width: 1180px;
        margin: 14px auto 0 auto;
      }
      .overview-spatial-card {
        padding: 14px 16px 10px 16px;
      }
      .overview-plot-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 12px;
        align-items: start;
      }
      .overview-plot-panel {
        min-width: 0;
      }
      .overview-plot-panel .shiny-plot-output {
        width: 100% !important;
      }
	      .overview-summary-grid {
	        display: grid;
	        grid-template-columns: repeat(3, minmax(0, 1fr));
	        gap: 14px;
	        align-items: stretch;
	      }
	      .overview-summary-grid .app-card:nth-child(1) {
	        border-top: 3px solid var(--app-accent);
	      }
	      .overview-summary-grid .app-card:nth-child(2) {
	        border-top: 3px solid var(--app-blue);
	      }
	      .overview-summary-grid .app-card:nth-child(3) {
	        border-top: 3px solid var(--app-violet);
	      }
      .overview-metric-row {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 10px;
        margin-bottom: 12px;
      }
	      .overview-metric {
	        background: var(--app-blue-soft);
	        border: 1px solid #d6e6ff;
	        border-radius: 8px;
	        padding: 11px 12px;
	      }
	      .overview-metric-value {
	        display: block;
	        color: var(--app-accent);
        font-size: 24px;
        font-weight: 700;
        line-height: 1.1;
      }
      .overview-metric-label {
        display: block;
        color: var(--app-muted);
        font-size: 12px;
        margin-top: 4px;
      }
      .overview-action-row {
        display: flex;
        justify-content: flex-end;
        margin-top: 14px;
      }
      .overview-action-row .btn {
        min-width: 170px;
      }
      .app-kv-item {
        display: flex;
        justify-content: space-between;
        gap: 12px;
        padding: 8px 0;
        border-bottom: 1px solid #eef3f8;
      }
      .app-kv-item:last-child {
        border-bottom: 0;
      }
      .app-kv-label {
        color: var(--app-muted);
      }
	      .app-kv-value {
	        color: var(--app-ink);
	        font-weight: 600;
	        text-align: right;
	      }
      .pca-shell {
        display: grid;
        grid-template-columns: minmax(280px, 340px) minmax(0, 1fr);
        gap: 18px;
        align-items: start;
        margin-top: 12px;
      }
      .pca-main {
        min-width: 0;
      }
      .pca-section {
        margin-bottom: 14px;
      }
      .pca-section > h4 {
        margin-top: 0;
        margin-bottom: 10px;
        font-size: 15px;
        font-weight: 700;
        color: #243447;
      }
      .pca-view-tabs .nav-tabs {
        border-bottom: 1px solid var(--app-border);
        margin-bottom: 14px;
      }
	      .pca-view-tabs .nav-tabs > li > a {
	        color: var(--app-muted);
	        font-weight: 700;
	        border-radius: 8px 8px 0 0;
	        border: 1px solid transparent;
	        padding: 9px 14px;
	      }
	      .pca-view-tabs .nav-tabs > li > a:hover,
	      .pca-view-tabs .nav-tabs > li > a:focus {
	        color: var(--app-accent);
	        background: var(--app-blue-soft);
	      }
	      .pca-view-tabs .nav-tabs > li.active > a,
	      .pca-view-tabs .nav-tabs > li.active > a:focus,
	      .pca-view-tabs .nav-tabs > li.active > a:hover {
	        color: var(--app-accent);
	        background: var(--app-accent-soft);
	        border-color: #d9e3ff #d9e3ff var(--app-accent-soft) #d9e3ff;
	      }
      .pca-view-tabs .tab-content {
        min-width: 0;
      }
      .pca-loading-controls {
        display: grid;
        grid-template-columns: repeat(3, minmax(170px, 1fr));
        gap: 12px;
        align-items: start;
        margin-bottom: 12px;
      }
      .pca-loading-controls .form-group {
        margin-bottom: 0;
      }
      .pca-orientation-help {
        background: var(--app-blue-soft);
        border: 1px solid #d6e6ff;
        border-radius: 8px;
        color: var(--app-ink);
        font-size: 12px;
        line-height: 1.45;
        margin-bottom: 10px;
        padding: 9px 10px;
      }
	      .pca-detail-grid {
	        display: grid;
	        grid-template-columns: minmax(360px, 0.95fr) minmax(0, 1.05fr);
	        gap: 14px;
	        align-items: stretch;
	      }
      .pca-detail-map-card .shiny-plot-output {
        width: 100% !important;
      }
      .pca-pc-row {
        padding: 9px 0;
        border-bottom: 1px solid #e8eef5;
      }
      .pca-pc-row:last-child {
        border-bottom: 0;
      }
      .pca-pc-row .checkbox {
        margin-top: 0;
        margin-bottom: 4px;
        font-weight: 700;
      }
      .pca-pc-row .form-group {
        margin-bottom: 0;
      }
	      .pca-pc-orientation .radio-inline {
	        margin-right: 10px;
	        color: var(--app-muted);
	      }
	      .pca-compact-summary {
	        padding: 5px 0;
	        border-bottom: 1px solid #eef3f8;
	      }
	      .pca-compact-summary:last-child {
	        border-bottom: 0;
	      }
	      .pca-summary-label {
	        color: var(--app-muted);
	        font-size: 12px;
	        margin-bottom: 5px;
	      }
	      .pca-selected-pills {
	        display: flex;
	        flex-wrap: nowrap;
	        gap: 4px;
	        align-items: center;
	        overflow: hidden;
	      }
		      .pca-pill {
	        display: inline-flex;
	        align-items: center;
	        justify-content: center;
		        border: 1px solid #c6d6f5;
		        background: var(--app-accent-soft);
		        border-radius: 8px;
		        color: var(--app-accent);
	        font-size: 12px;
	        font-weight: 700;
	        line-height: 1;
	        padding: 5px 7px;
	        white-space: nowrap;
	      }
	      .pca-summary-value {
	        color: #243447;
	        font-size: 12px;
	        font-weight: 600;
	        overflow-wrap: anywhere;
	      }
	      .pca-map-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
        gap: 14px;
      }
      .pca-map-grid.pc-count-1 {
        max-width: 560px;
        margin: 0 auto;
      }
      .pca-map-card {
        padding: 12px;
      }
      .pca-map-card .shiny-plot-output {
        width: 100% !important;
      }
      .pca-card-meta {
        color: var(--app-muted);
        font-size: 12px;
        margin-top: -5px;
        margin-bottom: 8px;
      }
      .pca-secondary-grid {
        display: grid;
        grid-template-columns: minmax(280px, 0.85fr) minmax(0, 1.15fr);
        gap: 14px;
      }
	      .loading-result-grid {
	        display: grid;
	        grid-template-columns: repeat(2, minmax(0, 1fr));
	        gap: 14px;
	        align-items: stretch;
	      }
      .loading-result-grid.single {
        grid-template-columns: minmax(280px, 680px);
        justify-content: center;
      }
	      .loading-panel {
	        background: #ffffff;
	        border: 1px solid #e5edf5;
	        border-top: 3px solid var(--app-blue);
	        border-radius: 8px;
	        padding: 12px;
	        min-width: 0;
	        height: 100%;
	      }
	      .loading-panel.negative {
	        border-top-color: var(--app-red);
	      }
	      .loading-panel h5 {
        margin-top: 0;
        margin-bottom: 10px;
	        color: var(--app-ink);
        font-weight: 700;
      }
      .loading-bar-row {
        padding: 7px 0;
        border-bottom: 1px solid #eef3f8;
      }
      .loading-bar-row:last-child {
        border-bottom: 0;
      }
      .loading-bar-head {
        display: flex;
        justify-content: space-between;
        gap: 12px;
        margin-bottom: 4px;
        font-size: 13px;
      }
	      .loading-gene {
	        color: var(--app-ink);
        font-weight: 600;
        overflow-wrap: anywhere;
      }
      .loading-value {
        color: var(--app-muted);
        font-variant-numeric: tabular-nums;
      }
	      .loading-bar-track {
	        height: 7px;
		        background: #edf2fb;
        border-radius: 999px;
        overflow: hidden;
      }
	      .loading-bar-fill {
	        height: 100%;
	        border-radius: 999px;
		        background: var(--app-blue);
	      }
	      .loading-bar-fill.negative {
	        background: var(--app-red);
	      }
      .mediation-shell {
        display: grid;
        grid-template-columns: minmax(280px, 340px) minmax(0, 1fr);
        gap: 18px;
        align-items: start;
        margin-top: 12px;
      }
      .mediation-main {
        min-width: 0;
      }
      .mediation-header {
        display: flex;
        justify-content: space-between;
        align-items: flex-start;
        gap: 14px;
        margin-bottom: 12px;
      }
	      .mediation-header h3 {
        margin: 0 0 4px 0;
        font-size: 20px;
        font-weight: 700;
	        color: var(--app-ink);
      }
      .mediation-header p {
        margin: 0;
      }
      .mediation-context {
        color: var(--app-muted);
        font-size: 13px;
        text-align: right;
      }
      .mediation-contrast-panel {
        margin-bottom: 14px;
      }
	      .mediation-contrast-panel h4 {
        margin-top: 0;
        margin-bottom: 10px;
        font-size: 15px;
        font-weight: 700;
	        color: var(--app-ink);
	      }
	      .effect-card-grid {
	        display: grid;
	        grid-template-columns: repeat(4, minmax(0, 1fr));
	        gap: 12px;
	        align-items: stretch;
	      }
	      .effect-card {
	        background: #ffffff;
	        border: 1px solid #e5edf5;
	        border-top: 3px solid var(--app-blue);
	        border-radius: 8px;
	        padding: 12px;
	        min-width: 0;
	        height: 100%;
	      }
	      .effect-card:nth-child(2) {
		        border-top-color: var(--app-violet);
	      }
	      .effect-card:nth-child(3) {
		        border-top-color: var(--app-accent);
	      }
	      .effect-card:nth-child(4) {
		        border-top-color: var(--app-slate-accent);
	      }
      .effect-label {
        color: var(--app-muted);
        font-size: 12px;
        margin-bottom: 6px;
      }
	      .effect-value {
	        color: var(--app-ink);
        font-size: 23px;
        font-weight: 700;
        line-height: 1.1;
        font-variant-numeric: tabular-nums;
      }
      .effect-interval,
      .effect-evidence {
        color: var(--app-muted);
        font-size: 12px;
        margin-top: 5px;
        font-variant-numeric: tabular-nums;
      }
      .mediation-result-grid {
        display: grid;
        grid-template-columns: minmax(0, 1fr);
        gap: 14px;
      }
      .mediation-view-tabs .nav-tabs {
        border-bottom: 1px solid var(--app-border);
        margin-bottom: 14px;
      }
	      .mediation-view-tabs .nav-tabs > li > a {
	        color: var(--app-muted);
	        font-weight: 700;
	        border-radius: 8px 8px 0 0;
	        border: 1px solid transparent;
	        padding: 9px 14px;
	      }
	      .mediation-view-tabs .nav-tabs > li > a:hover,
	      .mediation-view-tabs .nav-tabs > li > a:focus {
	        color: var(--app-accent);
	        background: var(--app-blue-soft);
	      }
	      .mediation-view-tabs .nav-tabs > li.active > a,
	      .mediation-view-tabs .nav-tabs > li.active > a:focus,
	      .mediation-view-tabs .nav-tabs > li.active > a:hover {
	        color: var(--app-accent);
	        background: var(--app-accent-soft);
	        border-color: #d9e3ff #d9e3ff var(--app-accent-soft) #d9e3ff;
	      }
      .mediation-view-tabs .tab-content {
        min-width: 0;
      }
      .mediation-path-table {
        width: 100%;
        border-collapse: collapse;
        font-size: 13px;
      }
      .mediation-path-table th {
        color: var(--app-muted);
        font-weight: 700;
        border-bottom: 1px solid #d9e2ec;
        padding: 8px 8px;
        text-align: left;
        white-space: nowrap;
      }
      .mediation-path-table td {
        border-bottom: 1px solid #eef3f8;
        padding: 8px 8px;
        vertical-align: top;
      }
      .mediation-path-table tr:last-child td {
        border-bottom: 0;
      }
      .mediation-header-actions {
        display: flex;
        justify-content: flex-end;
        align-items: center;
        gap: 8px;
        flex-shrink: 0;
      }
      .mediation-empty-action {
        margin-top: 14px;
      }
      .sensitivity-shell {
        display: grid;
        grid-template-columns: minmax(280px, 340px) minmax(0, 1fr);
        gap: 18px;
        align-items: start;
        margin-top: 12px;
      }
      .sensitivity-main {
        min-width: 0;
      }
      .sensitivity-header {
        display: flex;
        justify-content: space-between;
        align-items: flex-start;
        gap: 14px;
        margin-bottom: 12px;
      }
	      .sensitivity-header h3 {
        margin: 0 0 4px 0;
        font-size: 20px;
        font-weight: 700;
	        color: var(--app-ink);
      }
      .sensitivity-header p {
        margin: 0;
      }
      .sensitivity-context {
        color: var(--app-muted);
        font-size: 13px;
        text-align: right;
      }
      .sensitivity-result-grid {
        display: grid;
        grid-template-columns: minmax(0, 1fr);
        gap: 14px;
      }
      .sensitivity-view-tabs .nav-tabs {
        border-bottom: 1px solid var(--app-border);
        margin-bottom: 14px;
      }
	      .sensitivity-view-tabs .nav-tabs > li > a {
	        color: var(--app-muted);
	        font-weight: 700;
	        border-radius: 8px 8px 0 0;
	        border: 1px solid transparent;
	        padding: 9px 14px;
	      }
	      .sensitivity-view-tabs .nav-tabs > li > a:hover,
	      .sensitivity-view-tabs .nav-tabs > li > a:focus {
	        color: var(--app-accent);
	        background: var(--app-blue-soft);
	      }
	      .sensitivity-view-tabs .nav-tabs > li.active > a,
	      .sensitivity-view-tabs .nav-tabs > li.active > a:focus,
	      .sensitivity-view-tabs .nav-tabs > li.active > a:hover {
	        color: var(--app-accent);
	        background: var(--app-accent-soft);
	        border-color: #d9e3ff #d9e3ff var(--app-accent-soft) #d9e3ff;
	      }
      .sensitivity-view-tabs .tab-content {
        min-width: 0;
      }
      .sensitivity-tipping-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
        gap: 10px;
        margin-bottom: 12px;
      }
	      .sensitivity-tipping-card {
	        background: var(--app-slate-soft);
	        border: 1px solid #dce5f0;
	        border-top: 3px solid var(--app-slate-accent);
	        border-radius: 8px;
	        padding: 11px 12px;
	      }
      .sensitivity-tipping-label {
        color: var(--app-muted);
        font-size: 12px;
        margin-bottom: 5px;
      }
	      .sensitivity-tipping-value {
	        color: var(--app-accent);
	        font-size: 20px;
        font-weight: 700;
        font-variant-numeric: tabular-nums;
      }
      .sensitivity-decomp-panel {
        margin-bottom: 14px;
      }
	      .sensitivity-decomp-panel h4 {
        margin-top: 0;
        margin-bottom: 10px;
        font-size: 15px;
        font-weight: 700;
	        color: var(--app-ink);
      }
	      .sensitivity-effect-grid {
	        display: grid;
	        grid-template-columns: repeat(4, minmax(0, 1fr));
	        gap: 12px;
	        align-items: stretch;
	        margin-bottom: 12px;
	      }
      .sensitivity-mini-table {
        width: 100%;
        border-collapse: collapse;
        font-size: 13px;
      }
      .sensitivity-mini-table th {
        color: var(--app-muted);
        border-bottom: 1px solid #d9e2ec;
        padding: 7px 8px;
        text-align: left;
      }
      .sensitivity-mini-table td {
        border-bottom: 1px solid #eef3f8;
        padding: 7px 8px;
      }
      .sensitivity-mini-table tr:last-child td {
        border-bottom: 0;
      }
	      .sensitivity-warning {
	        color: #8a5a00;
	      }
      @media (max-width: 900px) {
        .setup-shell {
          grid-template-columns: 1fr;
        }
        .setup-sidebar {
          position: static;
        }
        .setup-main-grid {
          grid-template-columns: 1fr;
        }
	        .overview-shell {
	          width: 100%;
	        }
	        .tutorial-workflow,
	        .tutorial-tile-grid {
	          grid-template-columns: 1fr;
	        }
	        .tutorial-hero {
	          padding: 18px;
	        }
	        .tutorial-hero h2 {
	          font-size: 23px;
	        }
	        .overview-plot-grid,
        .overview-summary-grid,
        .overview-metric-row {
          grid-template-columns: 1fr;
        }
        .overview-action-row {
          justify-content: stretch;
        }
        .overview-action-row .btn {
          width: 100%;
        }
        .pca-shell,
        .pca-secondary-grid,
        .pca-loading-controls,
        .pca-detail-grid,
        .loading-result-grid,
        .loading-result-grid.single {
          grid-template-columns: 1fr;
        }
        .pca-map-grid.pc-count-1 {
          max-width: none;
        }
        .mediation-shell,
        .effect-card-grid {
          grid-template-columns: 1fr;
        }
        .mediation-header {
          display: block;
        }
        .mediation-header-actions {
          justify-content: flex-start;
          margin-top: 8px;
        }
        .mediation-context {
          text-align: left;
          margin-top: 6px;
        }
        .sensitivity-shell,
        .sensitivity-effect-grid {
          grid-template-columns: 1fr;
        }
        .sensitivity-header {
          display: block;
        }
        .sensitivity-context {
          text-align: left;
          margin-top: 6px;
        }
      }
      "))
	    )
	  ),
	  tabPanel(
	    "Tutorial",
	    div(
	      class = "tutorial-shell",
	      tags$section(
	        class = "tutorial-hero",
	        tags$div(class = "tutorial-eyebrow", "Guided workflow"),
	        tags$h2("Spatial transcriptomics mediation analysis, from gene expression to sensitivity checks"),
	        tags$p("Designed for spatial transcriptomics studies with a gene expression matrix, this app uses PCA-derived mediators to estimate mediation effects and examine how robust indirect effects are to mediator-outcome confounding.")
	      ),
	      tags$div(
	        class = "tutorial-workflow",
	        app_tutorial_step("1", "Data Setup", "Use the example data or upload matched spatial metadata and gene expression matrix files."),
	        app_tutorial_step("2", "Overview", "Check exposure groups, outcome patterns, model variables, and spatial distributions."),
	        app_tutorial_step("3", "PCA", "Choose which principal components will be used as mediators."),
	        app_tutorial_step("4", "Mediation", "Estimate direct, indirect, total, and mediator-specific effects."),
	        app_tutorial_step("5", "Sensitivity", "Assess how unmeasured mediator-outcome confounding could change the indirect effect.")
	      ),
	      tags$div(
	        class = "tutorial-tile-grid",
	        app_tutorial_tile(
	          "Default analysis",
	          "Bayesian mediation is the default mode for the spatial transcriptomics workflow. Frequentist mediation remains available as a comparison."
	        ),
	        app_tutorial_tile(
	          "What rho means",
	          "rho represents residual mediator-outcome correlation used in sensitivity analysis. rho = 0 matches the fitted mediation model."
	        ),
	        app_tutorial_tile(
	          "What to watch",
	          "Zero crossings show where an indirect effect reaches zero under a sensitivity scenario."
	        )
	      ),
	      app_tutorial_section(
	        "Data Requirements",
		        tags$p("The main dataset should contain one row per spatial observation, an observation ID, exposure, outcome, spatial coordinates, and optional covariates."),
		        tags$p("The mediator matrix is typically a gene expression matrix with the same observation ID and numeric gene or feature columns. Row order does not need to match because observations are matched by ID."),
	        open = TRUE
	      ),
	      app_tutorial_section(
	        "How the Workflow Fits Together",
	        tags$ol(
	          tags$li("Data Setup prepares the analysis object and matches observations across files."),
	          tags$li("Overview confirms that the prepared data and model variables look sensible."),
	          tags$li("PCA turns high-dimensional features into candidate mediator components."),
	          tags$li("Mediation decomposes the exposure-outcome association into direct and indirect components."),
	          tags$li("Sensitivity examines whether the indirect effect changes under residual mediator-outcome confounding.")
	        )
	      ),
	      app_tutorial_section(
	        "Sensitivity Modes",
	        tags$p("Common correlation applies the same residual correlation to all mediators at once."),
	        tags$p("One-mediator sensitivity varies one mediator at a time. The all-mediators comparison plot shows which mediator-specific indirect effects are most sensitive.")
	      )
	    )
	  ),
	  tabPanel(
	    "Data Setup",
    div(
      class = "setup-shell",
      tags$aside(
        class = "setup-sidebar",
        app_setup_section(
          "Data Source",
          radioButtons(
            "data_source",
            NULL,
            choices = c("Example Dataset" = "example", "Upload My Data" = "upload"),
            selected = "example"
          ),
          conditionalPanel(
            condition = "input.data_source == 'upload'",
            tags$div(
              class = "setup-meta",
              "Supported formats: CSV, TSV, TXT, XLSX, XLS."
            ),
            tags$hr(),
            h4("Main Analysis Data"),
            fileInput("main_file", "Upload main dataset", accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
            uiOutput("main_sheet_ui"),
            h4("Mediator Feature Matrix"),
            fileInput("feature_file", "Upload mediator feature matrix", accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
            uiOutput("feature_sheet_ui")
          )
        ),
        app_setup_section(
          "Variable Roles",
          uiOutput("observation_id_ui"),
          uiOutput("exposure_ui"),
          uiOutput("outcome_ui"),
          uiOutput("spatial_x_ui"),
          uiOutput("spatial_y_ui"),
          uiOutput("covariates_ui")
        ),
        app_setup_section(
          "Exposure Settings",
          uiOutput("reference_ui"),
          uiOutput("exposure_levels_ui")
        ),
        app_setup_section(
          "Observation Matching",
          uiOutput("feature_id_ui")
        ),
        app_setup_section(
          "Prepare",
          actionButton("prepare_analysis", "Prepare Analysis and Continue", class = "btn-primary"),
          uiOutput("prepare_message")
        )
	      ),
	      tags$main(
	        class = "setup-main",
	        uiOutput("data_setup_main_content")
	      )
	    )
  ),
	  tabPanel(
	    "Overview",
	    div(
	      class = "overview-shell",
	      app_panel_card(
	        "Spatial Overview",
        class = "overview-spatial-card",
        div(
          class = "overview-plot-grid",
          div(class = "overview-plot-panel", plotOutput("overview_exposure_spatial", height = 410)),
          div(class = "overview-plot-panel", plotOutput("overview_outcome_spatial", height = 410))
        )
      ),
      div(
        class = "overview-summary-grid",
        app_panel_card("Dataset", uiOutput("overview_dataset")),
        app_panel_card("Exposure Categories", uiOutput("overview_exposure")),
        app_panel_card("Model Setup", uiOutput("overview_variables"))
      ),
      div(
        class = "overview-action-row",
        actionButton("go_to_pca_from_overview", "Continue to PCA", class = "btn-primary")
      )
    )
  ),
  tabPanel(
    "PCA",
    div(
      class = "pca-shell",
      tags$aside(
        class = "setup-sidebar",
        app_setup_section(
          "PCs for Mediation",
          uiOutput("pca_pc_controls"),
          uiOutput("pca_selection_notice")
        ),
        app_setup_section(
          "Mediation Configuration",
          uiOutput("pca_confirmation_controls")
        )
      ),
	      tags$main(
	        class = "pca-main",
	        div(
          class = "pca-view-tabs",
          tabsetPanel(
            id = "pca_view",
            type = "tabs",
            tabPanel(
              "PC Detail",
              div(
                class = "pca-loading-controls",
                uiOutput("loading_pc_control"),
                selectInput("loading_n", "Number of loading genes/features", choices = c(10, 20, 30, 50), selected = 30),
                radioButtons(
                  "loading_sign",
                  "Show",
                  choices = c("Both" = "both", "Positive only" = "positive", "Negative only" = "negative"),
                  selected = "both"
                )
              ),
              div(
                class = "pca-detail-grid",
                app_panel_card(
                  "PC Score Map",
                  class = "pca-detail-map-card",
                  uiOutput("pca_detail_meta"),
                  plotOutput("pca_detail_spatial", height = 430)
                ),
                app_panel_card(
                  "Loadings",
                  uiOutput("loading_results")
                )
              )
            ),
            tabPanel(
              "Explained Variance",
              app_panel_card("Explained Variance", plotOutput("scree_plot", height = 430))
            )
          )
        )
      )
    )
	  ),
	  tabPanel(
	    "Mediation",
	    div(
	      class = "mediation-shell",
	      tags$aside(
	        class = "setup-sidebar",
	        uiOutput("mediation_input_context"),
	        app_setup_section(
	          "Mediation Settings",
	          radioButtons(
	            "analysis_framework",
	            "Framework",
	            choices = c("Bayesian" = "bayesian", "Frequentist" = "frequentist"),
	            selected = "bayesian"
	          ),
		          uiOutput("mediation_result_filters")
		        ),
	        conditionalPanel(
	          condition = "input.analysis_framework == 'frequentist'",
	          app_setup_section(
		            "Bootstrap CI",
		            numericInput("bootstrap_B", "Bootstrap replicates", value = 1000, min = 100, max = 5000, step = 100),
		            uiOutput("bootstrap_status")
		          )
		        )
	      ),
	      tags$main(
	        class = "mediation-main",
	        uiOutput("mediation_main_content")
	      )
	    )
	  ),
	  tabPanel(
	    "Sensitivity",
	    div(
	      class = "sensitivity-shell",
	      tags$aside(
	        class = "setup-sidebar",
	        uiOutput("sensitivity_baseline_context"),
	        app_setup_section(
	          "Sensitivity Settings",
	          radioButtons(
	            "sensitivity_framework",
	            "Framework",
	            choices = c("Bayesian" = "bayesian", "Frequentist" = "frequentist"),
	            selected = "bayesian"
	          ),
	          radioButtons(
	            "sensitivity_scenario",
	            "Sensitivity Type",
	            choices = c("Common correlation" = "common", "One mediator" = "one_at_a_time"),
	            selected = "common"
	          ),
	          uiOutput("sensitivity_contrast_control"),
	          uiOutput("sensitivity_mediator_control"),
	          uiOutput("sensitivity_rho_control")
	        )
	      ),
	      tags$main(
	        class = "sensitivity-main",
	        uiOutput("sensitivity_main_content")
	      )
	    )
	  )
)

server <- function(input, output, session) {
  analysis_state <- reactiveVal(NULL)
  setup_status <- reactiveVal(list(
    ok = FALSE,
    error = FALSE,
    message = "Click Prepare Analysis and Continue to start with the example dataset."
  ))
  bayesian_fit_state <- reactiveVal(NULL)
  bayesian_fit_cache <- reactiveVal(list())
  bayesian_fit_running <- reactiveVal(FALSE)
	  frequentist_fit_state <- reactiveVal(NULL)
	  frequentist_fit_cache <- reactiveVal(list())
	  frequentist_bootstrap_cache <- reactiveVal(list())
	  frequentist_sensitivity_cache <- reactiveVal(list())
  frequentist_fit_running <- reactiveVal(FALSE)
  sensitivity_warmup_status <- reactiveVal(NULL)
  confirmed_pca_spec <- reactiveVal(NULL)

  clear_workflow_results <- function() {
    confirmed_pca_spec(NULL)
    bayesian_fit_state(NULL)
    bayesian_fit_cache(list())
    frequentist_fit_state(NULL)
    frequentist_fit_cache(list())
    frequentist_bootstrap_cache(list())
    frequentist_sensitivity_cache(list())
    sensitivity_warmup_status(NULL)
  }

  set_analysis_unprepared <- function(message = "Upload both datasets, confirm variable roles, then prepare the analysis.", notify = FALSE) {
    clear_workflow_results()
    analysis_state(NULL)
    setup_status(list(ok = FALSE, error = FALSE, message = message))
    updateTabsetPanel(session, "main_nav", selected = "Data Setup")
    if (isTRUE(notify)) {
      showNotification(message, type = "message")
    }
  }

  safe_read_table <- function(file, sheet) {
    if (is.null(file)) {
      return(list(data = NULL, error = "No file uploaded."))
    }
    tryCatch(
      list(data = read_uploaded_table(file$datapath, sheet = sheet), error = NULL),
      error = function(e) list(data = NULL, error = conditionMessage(e))
    )
  }

  output$main_sheet_ui <- renderUI({
    req(input$main_file)
    ext <- tolower(tools::file_ext(input$main_file$name))
    if (!ext %in% c("xlsx", "xls")) {
      return(NULL)
    }
    sheets <- tryCatch(available_excel_sheets(input$main_file$datapath), error = function(e) character(0))
    if (length(sheets) == 0) {
      return(helpText("Excel sheet list is unavailable."))
    }
    selectInput("main_sheet", "Main dataset sheet", choices = sheets, selected = sheets[[1]])
  })

  output$feature_sheet_ui <- renderUI({
    req(input$feature_file)
    ext <- tolower(tools::file_ext(input$feature_file$name))
    if (!ext %in% c("xlsx", "xls")) {
      return(NULL)
    }
    sheets <- tryCatch(available_excel_sheets(input$feature_file$datapath), error = function(e) character(0))
    if (length(sheets) == 0) {
      return(helpText("Excel sheet list is unavailable."))
    }
    selectInput("feature_sheet", "Feature matrix sheet", choices = sheets, selected = sheets[[1]])
  })

  main_input <- reactive({
    if (identical(input$data_source %||% "example", "example")) {
      return(list(data = example_inputs$main_data, error = NULL, defaults = example_inputs$defaults))
    }
    x <- safe_read_table(input$main_file, input$main_sheet)
    x$defaults <- list()
    x
  })

  feature_input <- reactive({
    if (identical(input$data_source %||% "example", "example")) {
      return(list(data = example_inputs$feature_data, error = NULL, defaults = example_inputs$defaults))
    }
    x <- safe_read_table(input$feature_file, input$feature_sheet)
    x$defaults <- list()
    x
  })

  observeEvent(input$data_source, {
    if (identical(input$data_source %||% "example", "example")) {
      set_analysis_unprepared("Click Prepare Analysis and Continue to start with the example dataset.")
    } else {
      set_analysis_unprepared()
    }
  }, ignoreInit = TRUE)

  observeEvent(input$main_file, {
    if (identical(input$data_source %||% "example", "upload")) {
      set_analysis_unprepared("Main dataset selected. Confirm roles and prepare the analysis.")
    }
  }, ignoreInit = TRUE)

  observeEvent(input$feature_file, {
    if (identical(input$data_source %||% "example", "upload")) {
      set_analysis_unprepared("Feature matrix selected. Confirm matching and prepare the analysis.")
    }
  }, ignoreInit = TRUE)

  main_data_current <- reactive({
    x <- main_input()
    validate(need(is.null(x$error), x$error))
    x$data
  })

  feature_data_current <- reactive({
    x <- feature_input()
    validate(need(is.null(x$error), x$error))
    x$data
  })

  main_data_for_setup_ui <- reactive({
    x <- main_input()
    if (!is.null(x$error)) NULL else x$data
  })

  feature_data_for_setup_ui <- reactive({
    x <- feature_input()
    if (!is.null(x$error)) NULL else x$data
  })

  current_setup_spec <- reactive({
    main_data <- main_data_for_setup_ui()
    feature_data <- feature_data_for_setup_ui()
    if (is.null(main_data) || is.null(feature_data)) {
      return(NULL)
    }
    if (is.null(input$exposure_var) || is.null(input$outcome_var) ||
        is.null(input$spatial_x_var) || is.null(input$spatial_y_var) ||
        is.null(input$reference_category)) {
      return(NULL)
    }

    roles <- list(
      observation_id = input$observation_id %||% "",
      exposure = input$exposure_var,
      outcome = input$outcome_var,
      spatial_x = input$spatial_x_var,
      spatial_y = input$spatial_y_var,
      covariates = input$covariates %||% character(0),
      reference = input$reference_category
    )

    tryCatch(
      build_analysis_spec(
        data_source = input$data_source %||% "example",
        main_data = main_data,
        feature_data = feature_data,
        roles = roles,
        feature_id = input$feature_observation_id %||% "",
        excluded_feature_columns = character(0)
      ),
      error = function(e) structure(list(message = conditionMessage(e)), class = "setup_spec_error")
    )
  })

  setup_spec_is_error <- function(x) {
    inherits(x, "setup_spec_error")
  }

  default_role <- function(name, fallback = "") {
    if (identical(input$data_source %||% "example", "example")) {
      example_inputs$defaults[[name]] %||% fallback
    } else {
      fallback
    }
  }

  none_choices <- function(cols) {
    c("None" = "", cols)
  }

  output$main_data_summary <- renderTable({
    data <- main_data_current()
    numeric_cols <- numeric_column_names(data)
    data.frame(
      item = c("observations", "columns", "numeric columns"),
      value = c(nrow(data), ncol(data), length(numeric_cols)),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$feature_data_summary <- renderTable({
    data <- feature_data_current()
    numeric_cols <- numeric_column_names(data)
    data.frame(
      item = c("observations", "columns", "numeric columns"),
      value = c(nrow(data), ncol(data), length(numeric_cols)),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$main_data_preview <- renderDT({
    datatable(
      table_preview(main_data_current(), n_rows = 6, n_cols = 6),
      rownames = FALSE,
      class = "compact stripe nowrap",
      options = list(pageLength = 6, dom = "t", scrollX = TRUE, autoWidth = TRUE)
    )
  }, server = FALSE)

  output$feature_data_preview <- renderDT({
    datatable(
      table_preview(feature_data_current(), n_rows = 6, n_cols = 6),
      rownames = FALSE,
      class = "compact stripe nowrap",
      options = list(pageLength = 6, dom = "t", scrollX = TRUE, autoWidth = TRUE)
    )
  }, server = FALSE)

  output$main_column_types <- renderDT({
    datatable(column_type_table(main_data_current()), rownames = FALSE, options = list(pageLength = 8, dom = "tip"))
  }, server = FALSE)

  output$feature_numeric_status <- renderText({
    data <- feature_data_current()
    numeric_cols <- numeric_column_names(data)
    paste0(length(numeric_cols), " numeric columns detected before excluding ID/non-feature columns.")
  })

  output$observation_id_ui <- renderUI({
    data <- main_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selectInput("observation_id", "Observation ID", choices = none_choices(cols), selected = default_role("observation_id"))
  })

  output$exposure_ui <- renderUI({
    data <- main_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selectInput("exposure_var", "Exposure (X)", choices = cols, selected = default_role("exposure", cols[[1]]))
  })

  output$outcome_ui <- renderUI({
    data <- main_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selectInput("outcome_var", "Outcome (Y)", choices = cols, selected = default_role("outcome", cols[[1]]))
  })

  output$spatial_x_ui <- renderUI({
    data <- main_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selectInput("spatial_x_var", "Spatial X Coordinate", choices = cols, selected = default_role("spatial_x", cols[[1]]))
  })

  output$spatial_y_ui <- renderUI({
    data <- main_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selectInput("spatial_y_var", "Spatial Y Coordinate", choices = cols, selected = default_role("spatial_y", cols[[1]]))
  })

  output$covariates_ui <- renderUI({
    data <- main_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selected <- default_role("covariates", character(0))
    selectizeInput("covariates", "Additional Covariates", choices = cols, selected = selected, multiple = TRUE)
  })

  output$feature_id_ui <- renderUI({
    data <- feature_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selectizeInput(
      "feature_observation_id",
      "Feature Matrix Observation ID",
      choices = none_choices(cols),
      selected = default_role("feature_observation_id"),
      options = list(maxOptions = 2000)
    )
  })

  exposure_levels <- reactive({
    req(input$exposure_var)
    data <- main_data_for_setup_ui()
    req(data)
    unique(as.character(data[[input$exposure_var]]))
  })

  output$reference_ui <- renderUI({
    levels <- exposure_levels()
    if (length(levels) == 0) {
      return(NULL)
    }
    selected <- default_role("reference", levels[[1]])
    if (!selected %in% levels) {
      selected <- levels[[1]]
    }
    selectInput("reference_category", "Reference Category", choices = levels, selected = selected)
  })

  output$exposure_levels_ui <- renderUI({
    req(input$exposure_var)
    req(length(exposure_levels()) > 0)
    tags$p(
      class = "setup-meta",
      tags$strong("Observed exposure categories: "),
      paste(exposure_levels(), collapse = ", ")
    )
  })

  observeEvent(input$prepare_analysis, {
    roles <- list(
      observation_id = input$observation_id %||% "",
      exposure = input$exposure_var,
      outcome = input$outcome_var,
      spatial_x = input$spatial_x_var,
      spatial_y = input$spatial_y_var,
      covariates = input$covariates %||% character(0),
      reference = input$reference_category
    )

    result <- tryCatch({
      withProgress(message = "Preparing analysis data", value = 0, {
        setProgress(0.08, detail = "Checking selected variables and uploaded files")
        spec <- build_analysis_spec(
          data_source = input$data_source %||% "example",
          main_data = main_data_current(),
          feature_data = feature_data_current(),
          roles = roles,
          feature_id = input$feature_observation_id %||% "",
          excluded_feature_columns = character(0)
        )
        setProgress(0.35, detail = "Matching observations across datasets")
        setProgress(0.55, detail = "Preparing PCA inputs from the feature matrix")
        analysis_obj <- compute_analysis_from_spec(spec)
        setProgress(0.95, detail = "Preparing app views")
        list(ok = TRUE, spec = spec, analysis = analysis_obj, message = "Analysis data prepared successfully.")
      })
    }, error = function(e) {
      list(ok = FALSE, message = conditionMessage(e))
    })

    if (isTRUE(result$ok)) {
      analysis_state(result$analysis)
      setup_status(list(ok = TRUE, error = FALSE, message = result$message))
	      confirmed_pca_spec(NULL)
	      bayesian_fit_state(NULL)
	      bayesian_fit_cache(list())
	      frequentist_fit_state(NULL)
	      frequentist_fit_cache(list())
	      frequentist_bootstrap_cache(list())
	      frequentist_sensitivity_cache(list())
	      sensitivity_warmup_status(NULL)
      updateTabsetPanel(session, "main_nav", selected = "Overview")
    } else {
      setup_status(list(ok = FALSE, error = TRUE, message = result$message))
    }
  }, ignoreInit = TRUE)

  output$prepare_message <- renderUI({
    if (identical(input$data_source %||% "example", "upload") &&
        (is.null(input$main_file) || is.null(input$feature_file))) {
      return(tags$div(class = "setup-status setup-meta", "Upload both datasets, then prepare the analysis."))
    }
    status <- setup_status()
    color <- if (isTRUE(status$ok)) {
      "#1b6e3c"
    } else if (isTRUE(status$error)) {
      "#9b1c1c"
    } else {
      "#5f6f7f"
    }
    tags$div(class = "setup-status", style = paste0("color:", color, ";"), status$message)
  })

  output$validation_status <- renderTable({
    validation_status_table(analysis()$spec)
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  display_setup_spec <- reactive({
    spec <- current_setup_spec()
    if (is.null(spec) && identical(input$data_source %||% "example", "example")) {
      return(example_spec)
    }
    spec
  })

  output$data_setup_main_content <- renderUI({
    source <- input$data_source %||% "example"
    main <- main_input()
    feature <- feature_input()

    if (identical(source, "upload")) {
      no_main <- is.null(input$main_file)
      no_feature <- is.null(input$feature_file)
      if (no_main || no_feature) {
        return(app_panel_card(
          "Upload your data to begin",
          class = "setup-empty",
          tags$p("Use the controls on the left to upload two matched files:"),
          tags$ol(
            class = "setup-requirement-list",
            tags$li(
              tags$strong("A main analysis dataset"),
              tags$ul(
                tags$li("one observation ID column shared with the feature matrix"),
                tags$li("one categorical exposure/treatment column"),
                tags$li("one numeric outcome column"),
                tags$li("two numeric spatial coordinate columns"),
                tags$li("optional covariate columns")
              )
            ),
            tags$li(
              tags$strong("A mediator gene or feature matrix"),
              tags$ul(
                tags$li("the same observation ID column"),
                tags$li("numeric mediator, gene, or feature columns")
              )
            )
          ),
          tags$p("After both files are available, previews and validation will appear here.")
        ))
      }
    }

    read_errors <- c(main$error, feature$error)
    read_errors <- read_errors[!is.na(read_errors) & nzchar(read_errors)]
    if (length(read_errors) > 0) {
      return(app_panel_card(
        "Data upload issue",
        class = "setup-empty",
        tags$p(paste(unique(read_errors), collapse = " ")),
        tags$p("Check the files and sheet selections on the left.")
      ))
    }

    tagList(
      div(
        class = "setup-main-grid",
        app_panel_card(
        "Main Dataset Preview",
        tableOutput("main_data_summary"),
        div(class = "setup-preview-table", DTOutput("main_data_preview")),
        tags$details(
          tags$summary("Column types"),
          DTOutput("main_column_types")
        )
        ),
        app_panel_card(
        "Mediator Feature Matrix Preview",
        tableOutput("feature_data_summary"),
        div(class = "setup-preview-table", DTOutput("feature_data_preview")),
        tags$p(class = "setup-meta", textOutput("feature_numeric_status", inline = TRUE))
      )
      ),
      div(
        class = "setup-main-grid",
        app_panel_card(
          "Validation Summary",
          uiOutput("validation_status_items")
        ),
        app_panel_card(
          "Observation Matching",
          uiOutput("matching_summary")
        )
      )
    )
  })

  output$validation_status_items <- renderUI({
    spec <- display_setup_spec()
    if (is.null(spec)) {
      return(tags$p(class = "setup-meta", "Validation will appear after data files and variable roles are available."))
    }
    if (setup_spec_is_error(spec)) {
      return(tags$ul(
        class = "setup-validation-list",
        tags$li(tags$span(class = "setup-warning", "!"), spec$message)
      ))
    }

    rows <- validation_status_table(spec)
    tags$ul(
      class = "setup-validation-list",
      lapply(seq_len(nrow(rows)), function(i) {
        tags$li(
          tags$span(class = "setup-ok", HTML("&#10003;")),
          tags$strong(rows$Check[[i]]),
          ": ",
          rows$Status[[i]]
        )
      })
    )
  })

  output$matching_summary <- renderUI({
    spec <- display_setup_spec()
    if (is.null(spec)) {
      return(tags$p(class = "setup-meta", "Matching summary will appear after both datasets are available."))
    }
    if (setup_spec_is_error(spec)) {
      return(tags$p(class = "setup-meta", "Resolve the validation issue before matching can be summarized."))
    }

    tagList(
      tableOutput("matching_summary_table"),
      if (!is.null(spec$matching$warning)) {
        tags$p(class = "setup-meta", tags$span(class = "setup-warning", "!"), spec$matching$warning)
      } else {
        tags$p(class = "setup-meta", "Observation IDs are being used for matching.")
      }
    )
  })

  output$matching_summary_table <- renderTable({
    spec <- display_setup_spec()
    validate(need(!is.null(spec) && !setup_spec_is_error(spec), ""))
    unmatched <- (spec$matching$main_n - spec$matching$matched_n) +
      (spec$matching$feature_n - spec$matching$matched_n)
    data.frame(
      item = c("Main dataset", "Feature matrix", "Matched", "Unmatched"),
      value = c(spec$matching$main_n, spec$matching$feature_n, spec$matching$matched_n, unmatched),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  analysis <- reactive({
    analysis_state()
  })

  current_pca_spec <- reactive({
    selected <- selected_pc_ids(input, analysis())
    signs <- pc_sign_vector(input)
    list(
      selected_pcs = selected,
      pc_signs = signs,
      orientation_summary = pca_orientation_summary(selected, signs)
    )
  })

  pca_specs_equal <- function(a, b) {
    if (is.null(a) || is.null(b)) {
      return(FALSE)
    }
    if (!identical(a$selected_pcs, b$selected_pcs)) {
      return(FALSE)
    }
    selected <- a$selected_pcs
    identical(
      unname(a$pc_signs[selected]),
      unname(b$pc_signs[selected])
    )
  }

  pca_config_confirmed <- reactive({
    !is.null(confirmed_pca_spec())
  })

  pca_has_unconfirmed_changes <- reactive({
    spec <- confirmed_pca_spec()
    !is.null(spec) && !pca_specs_equal(spec, current_pca_spec())
  })

  confirmed_selected_pc_ids <- reactive({
    spec <- confirmed_pca_spec()
    if (is.null(spec)) {
      return(character(0))
    }
    spec$selected_pcs
  })

  confirmed_pc_sign_vector <- reactive({
    spec <- confirmed_pca_spec()
    if (is.null(spec)) {
      return(pc_sign_vector(input))
    }
    spec$pc_signs
  })

  mediation_sign_vector_confirmed <- reactive({
    pc_signs <- confirmed_pc_sign_vector()
    c(
      PC1_R = pc_signs[["PC1"]],
      PC2_R = pc_signs[["PC2"]],
      PC3 = pc_signs[["PC3"]]
    )
  })

  confirmed_pca_context <- reactive({
    spec <- confirmed_pca_spec()
    if (is.null(spec)) {
      return(NULL)
    }
    a <- analysis()
    list(
      exposure = a$spec$original_labels$exposure,
      reference = a$spec$exposure$reference,
      outcome = a$spec$original_labels$outcome,
      coordinates = paste(c(a$spec$original_labels$spatial_x, a$spec$original_labels$spatial_y), collapse = ", "),
      covariates = if (length(a$spec$covariates) == 0) "None" else paste(unname(a$spec$covariates), collapse = ", "),
      mediators = pca_spec_label(spec),
      orientation = spec$orientation_summary
    )
  })

  sensitivity <- reactive({
    freq <- frequentist_current_result()
    validate(need(!is.null(freq), "Run Frequentist Mediation first for the current data and selected PCs."))
    key <- freq$key
    cache <- frequentist_sensitivity_cache()
    if (!is.null(cache[[key]])) {
      return(cache[[key]])
    }
    value <- compute_frequentist_dynamic_sensitivity(
      frequentist_result = freq$result,
      rho_grid = bayesian_sensitivity_rho_grid
    )
    cache[[key]] <- value
    frequentist_sensitivity_cache(cache)
    value
  })

  set_sensitivity_warmup_status <- function(framework, state, message = NULL) {
    status <- sensitivity_warmup_status() %||% list()
    status[[framework]] <- list(state = state, message = message)
    sensitivity_warmup_status(status)
  }

  warm_bayesian_sensitivity_cache <- function(result, rho = 0) {
    if (is.null(result) || is.null(result$artifact)) {
      return(invisible(FALSE))
    }
    set_sensitivity_warmup_status("bayesian", "running")
    ok <- tryCatch({
      bayesian_sensitivity_curve_cache(
        artifact = result$artifact,
        artifact_key = result$key,
        mode = "common"
      )
      bayesian_sensitivity_selected_cache(
        artifact = result$artifact,
        artifact_key = result$key,
        rho = rho,
        mode = "common"
      )
      TRUE
    }, error = function(e) {
      set_sensitivity_warmup_status("bayesian", "error", conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok)) {
      set_sensitivity_warmup_status("bayesian", "ready")
    }
    invisible(ok)
  }

  warm_frequentist_sensitivity_cache <- function(result) {
    if (is.null(result) || is.null(result$result)) {
      return(invisible(FALSE))
    }
    set_sensitivity_warmup_status("frequentist", "running")
    ok <- tryCatch({
      key <- result$key
      cache <- frequentist_sensitivity_cache()
      if (is.null(cache[[key]])) {
        cache[[key]] <- compute_frequentist_dynamic_sensitivity(
          frequentist_result = result$result,
          rho_grid = bayesian_sensitivity_rho_grid
        )
        frequentist_sensitivity_cache(cache)
      }
      TRUE
    }, error = function(e) {
      set_sensitivity_warmup_status("frequentist", "error", conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok)) {
      set_sensitivity_warmup_status("frequentist", "ready")
    }
    invisible(ok)
  }

	  output$sensitivity_contrast_control <- renderUI({
	    a <- analysis()
	    contrast_choices <- c("All contrasts" = "all", stats::setNames(unname(a$contrast_labels), unname(a$contrast_labels)))
    selected <- input$sensitivity_contrast_filter %||% "all"
    if (!selected %in% unname(contrast_choices)) {
      selected <- "all"
    }
    selectInput("sensitivity_contrast_filter", "Contrast", choices = contrast_choices, selected = selected)
  })

	  contrast_label_current <- function(exposure_names) {
    labels <- analysis()$contrast_labels
    out <- unname(labels[exposure_names])
    ifelse(is.na(out), exposure_names, out)
  }

  filter_sensitivity_contrast <- function(df) {
    contrast <- input$sensitivity_contrast_filter %||% "all"
    if (identical(contrast, "all") || !"contrast" %in% names(df)) {
      return(df)
    }
    dplyr::filter(df, contrast == !!contrast)
  }

  bayesian_artifact_available <- reactive({
    if (!isTRUE(pca_config_confirmed())) {
      return(FALSE)
    }
    isTRUE(analysis()$bayesian_artifact_allowed) &&
      validated_selection_selected(confirmed_selected_pc_ids())
  })

  bayesian_current_mediators <- reactive({
    if (isTRUE(bayesian_artifact_available())) {
      return(mediators)
    }
    confirmed_selected_pc_ids()
  })

  bayesian_dynamic_settings <- reactive({
    bayesian_dynamic_default_settings()
  })

  bayesian_current_cache_key <- reactive({
    req(pca_config_confirmed())
    bayesian_dynamic_cache_key(
      analysis(),
      selected_pcs = confirmed_selected_pc_ids(),
      pc_signs = confirmed_pc_sign_vector(),
      settings = bayesian_dynamic_settings()
    )
  })

  bayesian_current_result <- reactive({
    if (!isTRUE(pca_config_confirmed())) {
      return(NULL)
    }
    if (isTRUE(bayesian_artifact_available())) {
      return(list(source = "validated_artifact", key = "emory_validated", artifact = bayesian_app_artifact))
    }

    key <- bayesian_current_cache_key()
    state <- bayesian_fit_state()
    if (!is.null(state) && identical(state$key, key)) {
      return(state)
    }

    cached <- bayesian_fit_cache()[[key]]
    if (!is.null(cached)) {
      bayesian_fit_state(cached)
      return(cached)
    }

    NULL
  })

  bayesian_current_artifact <- reactive({
    result <- bayesian_current_result()
    if (is.null(result)) NULL else result$artifact
  })

  pc_signs_are_default <- reactive({
    if (!isTRUE(pca_config_confirmed())) {
      return(FALSE)
    }
    all(confirmed_pc_sign_vector()[confirmed_selected_pc_ids()] == 1)
  })

  frequentist_example_default_available <- reactive({
    if (!isTRUE(pca_config_confirmed())) {
      return(FALSE)
    }
    isTRUE(analysis()$bayesian_artifact_allowed) &&
      validated_selection_selected(confirmed_selected_pc_ids()) &&
      isTRUE(pc_signs_are_default())
  })

  frequentist_current_mediators <- reactive({
    if (isTRUE(frequentist_example_default_available())) {
      return(validated_pc_ids)
    }
    confirmed_selected_pc_ids()
  })

  frequentist_current_cache_key <- reactive({
    req(pca_config_confirmed())
    if (isTRUE(frequentist_example_default_available())) {
      return("emory_frequentist_default")
    }
    frequentist_dynamic_cache_key(
      analysis(),
      selected_pcs = confirmed_selected_pc_ids(),
      pc_signs = confirmed_pc_sign_vector()
    )
  })

  frequentist_current_result <- reactive({
    if (!isTRUE(pca_config_confirmed())) {
      return(NULL)
    }
    if (isTRUE(frequentist_example_default_available())) {
      return(list(
        source = "validated_example",
        key = frequentist_current_cache_key(),
        result = make_frequentist_example_fit(analysis())
      ))
    }

    key <- frequentist_current_cache_key()
    state <- frequentist_fit_state()
    if (!is.null(state) && identical(state$key, key)) {
      return(state)
    }

    cached <- frequentist_fit_cache()[[key]]
    if (!is.null(cached)) {
      frequentist_fit_state(cached)
      return(cached)
    }

    NULL
  })

	  frequentist_unavailable_message <- function() {
	    if (!isTRUE(pca_config_confirmed())) {
	      return("Confirm mediator configuration in PCA before running Frequentist Mediation.")
	    }
	    "Click Run Frequentist Mediation to fit the frequentist model for the confirmed prepared data, selected PCs, and sign orientations."
	  }

	  mediation_result_available <- reactive({
	    framework <- input$analysis_framework %||% "bayesian"
	    if (identical(framework, "bayesian")) {
	      return(!is.null(bayesian_current_artifact()))
	    }
	    !is.null(frequentist_current_result())
	  })

	  output$mediation_result_filters <- renderUI({
	    if (!isTRUE(mediation_result_available())) {
	      return(NULL)
	    }

	    a <- analysis()
	    contrast_choices <- c("All contrasts" = "all", stats::setNames(unname(a$contrast_labels), unname(a$contrast_labels)))
	    selected_contrast <- input$mediation_contrast_filter %||% "all"
	    if (!selected_contrast %in% unname(contrast_choices)) {
	      selected_contrast <- "all"
	    }

	    mediator_ids <- if (identical(input$analysis_framework %||% "bayesian", "frequentist")) {
	      frequentist_current_mediators()
	    } else {
	      bayesian_current_mediators()
	    }
	    mediator_choices <- mediator_filter_choices(mediator_ids)
	    selected_mediator <- input$mediation_mediator_filter %||% "all"
	    if (!selected_mediator %in% unname(mediator_choices)) {
	      selected_mediator <- "all"
	    }

	    tagList(
	      selectInput("mediation_contrast_filter", "Contrast", choices = contrast_choices, selected = selected_contrast),
	      selectInput("mediation_mediator_filter", "Mediator", choices = mediator_choices, selected = selected_mediator)
	    )
	  })

	  sensitivity_mediator_ids <- reactive({
    if (identical(input$sensitivity_framework %||% "bayesian", "bayesian")) {
      artifact <- bayesian_current_artifact()
      if (!is.null(artifact)) {
        return(bayes_A_bundle_mediators(artifact$sensitivity_A$residual_draws))
      }
      return(bayesian_current_mediators())
    }
    frequentist_current_mediators()
  })

  output$sensitivity_mediator_control <- renderUI({
    if (!identical(input$sensitivity_scenario %||% "common", "one_at_a_time")) {
      return(NULL)
    }
    choices <- mediator_filter_choices(sensitivity_mediator_ids())[-1]
    if (length(choices) == 0) {
      return(tags$p(class = "setup-meta", "Confirm mediators and fit mediation first."))
    }
    selected <- input$sensitivity_one_mediator %||% unname(choices)[[1]]
    if (!selected %in% unname(choices)) {
      selected <- unname(choices)[[1]]
    }
    selectInput("sensitivity_one_mediator", "Mediator", choices = choices, selected = selected)
  })

  output$sensitivity_rho_control <- renderUI({
    sliderInput(
      "sensitivity_active_rho",
      "Residual correlation (rho)",
      min = min(bayesian_sensitivity_rho_grid),
      max = max(bayesian_sensitivity_rho_grid),
      value = input$sensitivity_active_rho %||% 0,
      step = 0.05
    )
  })

  bayesian_unavailable_message <- function() {
    if (!isTRUE(pca_config_confirmed())) {
      return("Confirm mediator configuration in PCA before running Bayesian Mediation.")
    }
    "Click Run Bayesian Mediation to fit the Bayesian model for the confirmed prepared data, selected PCs, and sign orientations."
  }

  output$bayesian_mediation_guard <- renderUI({
    if (isTRUE(bayesian_artifact_available())) {
      return(tags$div(style = "color:#1b6e3c;", "Bayesian results are available for the current selection."))
    }
    tagList(
      fluidRow(
        column(
          width = 4,
          actionButton(
            "run_bayesian_mediation",
            "Run Bayesian Mediation",
            class = "btn-primary"
          )
        ),
        column(width = 8, uiOutput("bayesian_fit_status"))
      ),
      tags$div(style = "color:#8a5a00;", bayesian_unavailable_message())
    )
  })

  output$bayesian_sensitivity_guard <- renderUI({
    if (!is.null(bayesian_current_artifact())) {
      return(NULL)
    }
    tags$div(style = "color:#9b1c1c;", tags$strong("Bayesian sensitivity unavailable: "), "Run Bayesian Mediation first to enable sensitivity analysis.")
  })

  output$frequentist_sensitivity_guard <- renderUI({
    if (!is.null(frequentist_current_result())) {
      return(NULL)
    }
    tags$div(style = "color:#9b1c1c;", tags$strong("Frequentist sensitivity unavailable: "), "Run Frequentist Mediation first for the current data and selected PCs.")
  })

  output$mediation_input_context <- renderUI({
    context <- confirmed_pca_context()
    if (is.null(context)) {
      return(tags$div(
        class = "workflow-context-card",
        h5("Analysis Input"),
        tags$p(class = "workflow-note", "Confirm mediators in PCA before running mediation.")
      ))
    }

    tagList(
      tags$div(
        class = "workflow-context-card",
        h5("Analysis Input"),
        app_kv_item("Exposure", context$exposure),
        app_kv_item("Reference", context$reference),
        app_kv_item("Mediators", context$mediators),
        app_kv_item("Covariates", context$covariates)
      ),
      if (isTRUE(pca_has_unconfirmed_changes())) {
        tags$p(class = "workflow-note warning", "PCA edits are unconfirmed. Mediation still uses the confirmed mediator configuration.")
      }
    )
  })

  output$sensitivity_baseline_context <- renderUI({
    framework <- input$sensitivity_framework %||% "bayesian"
    context <- confirmed_pca_context()
    contrast <- input$sensitivity_contrast_filter %||% "all"
    contrast_label_ui <- if (identical(contrast, "all")) "All contrasts" else contrast

    if (is.null(context)) {
      return(tags$div(
        class = "workflow-context-card",
        h5("Analysis Input"),
        tags$p(class = "workflow-note", "Confirm mediators in PCA before running sensitivity analysis.")
      ))
    }

    tags$div(
      class = "workflow-context-card",
      h5("Analysis Input"),
      app_kv_item("Framework", if (identical(framework, "bayesian")) "Bayesian" else "Frequentist"),
      app_kv_item("Contrast", contrast_label_ui),
      app_kv_item("Mediators", context$mediators)
    )
  })

	  observeEvent(input$continue_to_sensitivity, {
	    updateRadioButtons(
	      session,
	      "sensitivity_framework",
	      selected = input$analysis_framework %||% "bayesian"
	    )
	    updateTabsetPanel(session, "main_nav", selected = "Sensitivity")
	  }, ignoreInit = TRUE)

  observeEvent(input$go_to_mediation_from_sensitivity, {
    updateTabsetPanel(session, "main_nav", selected = "Mediation")
  }, ignoreInit = TRUE)

  observeEvent(input$go_to_pca_from_overview, {
    updateTabsetPanel(session, "main_nav", selected = "PCA")
  }, ignoreInit = TRUE)

  observe({
	    data_done <- isTRUE(setup_status()$ok) && !is.null(analysis())
	    pca_done <- isTRUE(pca_config_confirmed()) && !isTRUE(pca_has_unconfirmed_changes())
	    mediation_done <- pca_done && (
	      !is.null(bayesian_current_artifact()) ||
	        !is.null(frequentist_current_result())
	    )
	    sensitivity_available <- mediation_done

    steps <- list(
      list(value = "Data Setup", workflow = TRUE, completed = data_done, available = TRUE),
      list(value = "Overview", workflow = FALSE, completed = data_done, available = data_done),
      list(value = "PCA", workflow = TRUE, completed = pca_done, available = data_done),
      list(value = "Mediation", workflow = TRUE, completed = mediation_done, available = pca_done),
      list(value = "Sensitivity", workflow = TRUE, completed = sensitivity_available, available = sensitivity_available)
    )
    session$sendCustomMessage("setNavWorkflow", list(steps = steps))
  })

  output$bayesian_fit_status <- renderUI({
    if (isTRUE(bayesian_fit_running())) {
      return(tags$div(style = "color:#8a5a00;", "Bayesian model fitting is running."))
    }
    if (!is.null(bayesian_current_artifact())) {
      result <- bayesian_current_result()
      diagnostics <- result$artifact$diagnostics %||% NULL
      status <- if (!is.null(diagnostics) && all(diagnostics$overall_pass, na.rm = TRUE)) {
        "Bayesian fit is complete and diagnostics passed."
      } else if (!is.null(diagnostics)) {
        "Bayesian fit is complete with diagnostic warnings."
      } else {
        "Bayesian results are ready."
      }
      return(tags$div(style = "color:#1b6e3c;", status))
    }
    tags$div(style = "color:#666;", "No Bayesian fit has been run for the current setup.")
  })

  output$bayesian_diagnostics_notice <- renderUI({
    artifact <- bayesian_current_artifact()
    if (is.null(artifact) || is.null(artifact$diagnostics)) {
      return(NULL)
    }
    diagnostics <- artifact$diagnostics
    color <- if (all(diagnostics$overall_pass, na.rm = TRUE)) "#1b6e3c" else "#8a5a00"
    tags$details(
      tags$summary("Bayesian fit diagnostics"),
      tags$div(
        style = paste0("color:", color, "; margin: 6px 0;"),
        if (all(diagnostics$overall_pass, na.rm = TRUE)) {
          "All monitored MCMC diagnostics passed."
        } else {
          "One or more monitored MCMC diagnostics did not pass. Interpret posterior summaries with caution."
        }
      ),
      tableOutput("bayesian_diagnostics_table")
    )
  })

  output$bayesian_diagnostics_table <- renderTable({
    req(bayesian_current_artifact())
    req(bayesian_current_artifact()$diagnostics)
    bayesian_current_artifact()$diagnostics |>
      dplyr::transmute(
        model,
        max_Rhat = round(max_Rhat, 3),
        min_bulk_ESS = round(min_bulk_ESS, 0),
        min_tail_ESS = round(min_tail_ESS, 0),
        divergent_transitions,
        max_treedepth_hits,
        overall_pass
      )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  observeEvent(input$run_bayesian_mediation, {
    if (!isTRUE(pca_config_confirmed())) {
      showNotification("Confirm mediator configuration in PCA before running Bayesian mediation.", type = "warning")
      return(NULL)
    }
    req(!isTRUE(bayesian_artifact_available()))
    if (isTRUE(bayesian_fit_running())) {
      showNotification("Bayesian mediation is already running for this session.", type = "warning")
      return(NULL)
    }

    key <- bayesian_current_cache_key()
    cached <- bayesian_fit_cache()[[key]]
    if (!is.null(cached)) {
      bayesian_fit_state(cached)
      showNotification("Reused Bayesian mediation results for the current setup.", type = "message")
      warm_bayesian_sensitivity_cache(cached)
      return(NULL)
    }

    bayesian_fit_running(TRUE)
    session$sendCustomMessage("setDisabled", list(id = "run_bayesian_mediation", disabled = TRUE))
    on.exit({
      bayesian_fit_running(FALSE)
      session$sendCustomMessage("setDisabled", list(id = "run_bayesian_mediation", disabled = FALSE))
    }, add = TRUE)

    settings <- bayesian_dynamic_settings()
    result <- tryCatch({
      artifact <- withProgress(message = "Running Bayesian mediation model", value = 0, {
        setProgress(0.02, detail = "Preparing model data")
        progress_values <- c(
          prepare = 0.08,
          mediator_model = 0.20,
          outcome_model = 0.58,
          posterior_processing = 0.84,
          diagnostics = 0.94
        )
        progress_labels <- c(
          prepare = "Preparing model data",
          mediator_model = "Fitting mediator model (stage 1 of 2)",
          outcome_model = "Fitting outcome model (stage 2 of 2)",
          posterior_processing = "Summarizing posterior draws",
          diagnostics = "Checking MCMC diagnostics"
        )
        make_bayesian_dynamic_artifact(
          analysis_obj = analysis(),
          selected_pcs = confirmed_selected_pc_ids(),
          pc_signs = confirmed_pc_sign_vector(),
          settings = settings,
          progress_callback = function(step, detail = NULL) {
            stage_label <- if (step %in% names(progress_labels)) progress_labels[[step]] else step
            if (!is.null(detail) && nzchar(detail)) {
              stage_label <- paste(stage_label, detail, sep = ": ")
            }
            stage_value <- if (step %in% names(progress_values)) progress_values[[step]] else 0.5
            setProgress(stage_value, detail = stage_label)
          }
        )
      })
      list(ok = TRUE, artifact = artifact)
    }, error = function(e) {
      list(ok = FALSE, message = conditionMessage(e))
    })

    if (isTRUE(result$ok)) {
      state <- list(source = "dynamic_fit", key = key, artifact = result$artifact)
      cache <- bayesian_fit_cache()
      cache[[key]] <- state
      bayesian_fit_cache(cache)
      bayesian_fit_state(state)
      showNotification("Bayesian mediation fit completed.", type = "message")
      warm_bayesian_sensitivity_cache(state)
    } else {
      showNotification(paste("Bayesian mediation failed:", result$message), type = "error", duration = 12)
    }
  }, ignoreInit = TRUE)

  output$frequentist_mediation_guard <- renderUI({
    if (isTRUE(frequentist_example_default_available())) {
      return(tags$div(style = "color:#1b6e3c;", "Frequentist results are available for the current selection."))
    }
    tagList(
      fluidRow(
        column(
          width = 4,
          actionButton(
            "run_frequentist_mediation",
            "Run Frequentist Mediation",
            class = "btn-primary"
          )
        ),
        column(width = 8, uiOutput("frequentist_fit_status"))
      ),
      tags$div(style = "color:#8a5a00;", frequentist_unavailable_message())
    )
  })

  output$frequentist_fit_status <- renderUI({
    if (isTRUE(frequentist_fit_running())) {
      return(tags$div(style = "color:#8a5a00;", "Frequentist mediation is running."))
    }
    if (!is.null(frequentist_current_result())) {
      return(tags$div(style = "color:#1b6e3c;", "Frequentist mediation fit is ready."))
    }
    tags$div(style = "color:#666;", "No frequentist fit has been run for the current setup.")
  })

	  observeEvent(input$run_frequentist_mediation, {
    if (!isTRUE(pca_config_confirmed())) {
      showNotification("Confirm mediator configuration in PCA before running Frequentist mediation.", type = "warning")
      return(NULL)
    }
    if (isTRUE(frequentist_fit_running())) {
      showNotification("Frequentist mediation is already running for this session.", type = "warning")
      return(NULL)
    }

    key <- frequentist_current_cache_key()
    cached <- frequentist_fit_cache()[[key]]
    if (!is.null(cached)) {
      frequentist_fit_state(cached)
      showNotification("Reused frequentist mediation results for the current setup.", type = "message")
      warm_frequentist_sensitivity_cache(cached)
      return(NULL)
    }

    frequentist_fit_running(TRUE)
    session$sendCustomMessage("setDisabled", list(id = "run_frequentist_mediation", disabled = TRUE))
    on.exit({
      frequentist_fit_running(FALSE)
      session$sendCustomMessage("setDisabled", list(id = "run_frequentist_mediation", disabled = FALSE))
    }, add = TRUE)

    result <- tryCatch({
      fit <- make_frequentist_dynamic_fit(
        analysis_obj = analysis(),
        selected_pcs = confirmed_selected_pc_ids(),
        pc_signs = confirmed_pc_sign_vector()
      )
      list(ok = TRUE, fit = fit)
    }, error = function(e) {
      list(ok = FALSE, message = conditionMessage(e))
    })

    if (isTRUE(result$ok)) {
      state <- list(source = "dynamic_fit", key = key, result = result$fit)
      cache <- frequentist_fit_cache()
      cache[[key]] <- state
      frequentist_fit_cache(cache)
      frequentist_fit_state(state)
      showNotification("Frequentist mediation fit completed.", type = "message")
      warm_frequentist_sensitivity_cache(state)
    } else {
      showNotification(paste("Frequentist mediation failed:", result$message), type = "error", duration = 12)
    }
  }, ignoreInit = TRUE)

  overview_spatial_df <- reactive({
    analysis()$spatial_inputs$plot_df
  })

  overview_spatial_limits <- reactive({
    spatial_plot_limits(overview_spatial_df())
  })

  output$overview_exposure_spatial <- renderPlot({
    spec <- analysis()$spec
    limits <- overview_spatial_limits()
    overview_spatial_df() |>
      ggplot(aes(x = x_plot, y = y_plot, color = X_cat)) +
      geom_point(size = 1.35, alpha = 0.85) +
      coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
      labs(
        title = paste("Exposure Spatial Distribution:", spec$original_labels$exposure),
        x = spec$original_labels$spatial_x,
        y = spec$original_labels$spatial_y,
        color = spec$original_labels$exposure
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        plot.margin = margin(4, 4, 4, 4),
        legend.margin = margin(0, 0, 0, 0)
      )
  })

  output$overview_outcome_spatial <- renderPlot({
    spec <- analysis()$spec
    limits <- overview_spatial_limits()
    overview_spatial_df() |>
      ggplot(aes(x = x_plot, y = y_plot, color = Y)) +
      geom_point(size = 1.35, alpha = 0.9) +
      coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
      scale_color_gradient(low = "#2C7BB6", high = "#D7191C") +
      labs(
        title = paste("Outcome Spatial Distribution:", spec$original_labels$outcome),
        x = spec$original_labels$spatial_x,
        y = spec$original_labels$spatial_y,
        color = spec$original_labels$outcome
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        plot.margin = margin(4, 4, 4, 4),
        legend.margin = margin(0, 0, 0, 0)
      )
  })

  output$overview_dataset <- renderUI({
    a <- analysis()
    tags$div(
      class = "overview-metric-row",
      tags$div(
        class = "overview-metric",
        tags$span(class = "overview-metric-value", format(nrow(a$dat_full_allpc), big.mark = ",")),
        tags$span(class = "overview-metric-label", "observations")
      ),
      tags$div(
        class = "overview-metric",
        tags$span(class = "overview-metric-value", format(a$spec$mediator_features$n_features, big.mark = ",")),
        tags$span(class = "overview-metric-label", "mediator features")
      )
    )
  })

  output$overview_exposure <- renderUI({
    counts <- analysis()$dat_full_allpc |>
      count(X_cat, name = "spots") |>
      rename(category = X_cat)

    tags$div(
      lapply(seq_len(nrow(counts)), function(i) {
        app_kv_item(
          as.character(counts$category[[i]]),
          format(counts$spots[[i]], big.mark = ",")
        )
      })
    )
  })

  output$overview_variables <- renderUI({
    a <- analysis()
    tags$div(
      app_kv_item("Exposure", a$spec$original_labels$exposure),
      app_kv_item("Reference", a$spec$exposure$reference),
      app_kv_item("Outcome", a$spec$original_labels$outcome),
      app_kv_item("Coordinates", paste(c(a$spec$original_labels$spatial_x, a$spec$original_labels$spatial_y), collapse = ", ")),
      app_kv_item(
        "Additional Covariates",
        if (length(a$spec$covariates) == 0) "none" else paste(unname(a$spec$covariates), collapse = ", ")
      )
    )
  })

  output$scree_plot <- renderPlot({
    analysis()$scree_df |>
      slice_head(n = 30) |>
      ggplot(aes(x = PC, y = var_explained)) +
      geom_col(fill = "#4C78A8") +
      geom_line(aes(y = cum_var_explained), color = "#F58518", linewidth = 0.8) +
      geom_point(aes(y = cum_var_explained), color = "#F58518", size = 1.6) +
      labs(x = "Principal component", y = "Proportion variance explained") +
      theme_bw(base_size = 12)
  })

  observeEvent(input$confirm_pca_mediators, {
    spec <- current_pca_spec()
    if (length(spec$selected_pcs) == 0) {
      showNotification("Select at least one PC before confirming mediators.", type = "warning")
      return(NULL)
    }

    previous <- confirmed_pca_spec()
    changed <- !pca_specs_equal(previous, spec)
    confirmed_pca_spec(spec)
    if (isTRUE(changed)) {
	      bayesian_fit_state(NULL)
	      frequentist_fit_state(NULL)
	      frequentist_bootstrap_cache(list())
	      sensitivity_warmup_status(NULL)
    }
    bayes_ready <- bayesian_current_result()
    if (!is.null(bayes_ready)) {
      warm_bayesian_sensitivity_cache(bayes_ready)
    }
    freq_ready <- frequentist_current_result()
    if (!is.null(freq_ready)) {
      warm_frequentist_sensitivity_cache(freq_ready)
    }
    showNotification("Mediator configuration confirmed.", type = "message")
    updateTabsetPanel(session, "main_nav", selected = "Mediation")
  }, ignoreInit = TRUE)

	  output$pca_confirmation_controls <- renderUI({
	    spec <- current_pca_spec()
	    confirmed <- confirmed_pca_spec()

    status <- if (is.null(confirmed)) {
      tags$p(class = "workflow-note", "Confirm this configuration before running mediation.")
    } else if (isTRUE(pca_has_unconfirmed_changes())) {
      tags$p(class = "workflow-note warning", "PC selection has changed. Confirm again to update the mediation analysis.")
    } else {
      tags$p(class = "workflow-note ok", HTML("&#10003; Mediator configuration confirmed"))
    }

    tagList(
	      tags$div(
	        class = "workflow-context-card",
	        h5("Current PCA selection"),
	        tags$div(
	          class = "pca-compact-summary",
	          tags$div(class = "pca-summary-label", "Selected mediators"),
	          if (length(spec$selected_pcs) == 0) {
	            tags$div(class = "pca-summary-value", "None")
	          } else {
	            tags$div(
	              class = "pca-selected-pills",
	              lapply(spec$selected_pcs, function(pc) tags$span(class = "pca-pill", pc))
	            )
	          }
	        ),
	        tags$div(
	          class = "pca-compact-summary",
	          tags$div(class = "pca-summary-label", "Orientation"),
	          tags$div(class = "pca-summary-value", spec$orientation_summary)
	        )
	      ),
      status,
      tags$div(
        class = "workflow-action",
        actionButton("confirm_pca_mediators", "Confirm Mediators & Continue", class = "btn-primary")
      )
    )
  })

  output$pca_selection_notice <- renderUI({
    selected <- selected_pc_ids(input, analysis())
    if (length(selected) == 0) {
      return(tags$p(class = "setup-meta", tags$span(class = "setup-warning", "!"), "Select at least one PC to display maps and loading genes."))
    }
    if (isTRUE(bayesian_artifact_available())) {
      return(tags$p(class = "setup-meta", "Bayesian results are available for PC1, PC2, PC3 with the current orientation."))
    }

    tags$div(
      class = "setup-meta",
      style = "color:#8a5a00;",
      tags$strong("Fit required: "),
      "the selected PCs and sign orientations will be used when you run Bayesian or Frequentist mediation for the current setup."
    )
  })

  output$pca_pc_controls <- renderUI({
    candidates <- pc_candidate_ids(analysis(), max_pcs = 6)

    tagList(
      tags$div(
        class = "pca-orientation-help",
        tags$strong("Sign orientation: "),
        "Keep/Reverse only flips the sign direction of a PC. Variance explained is unchanged; the selected orientation is used for maps, loading direction, and mediation."
      ),
      lapply(seq_along(candidates), function(i) {
        pc <- candidates[[i]]
        selected_id <- paste0("pca_pc", i, "_selected")
        orientation_id <- paste0("pc", i, "_orientation")
        is_selected <- input[[selected_id]] %||% TRUE

        tags$div(
          class = "pca-pc-row",
          checkboxInput(selected_id, pc, value = isTRUE(is_selected)),
          if (isTRUE(is_selected)) {
            tags$div(
              class = "pca-pc-orientation",
              radioButtons(
                orientation_id,
                NULL,
                choices = c("Keep" = "keep", "Reverse" = "reverse"),
                selected = input[[orientation_id]] %||% "keep",
                inline = TRUE
              )
            )
          }
        )
      })
    )
  })

  oriented_spatial_df <- reactive({
    orient_spatial_scores(
      analysis()$spatial_inputs$plot_df,
      pc_sign_vector(input),
      pca_fit = analysis()$pca$pca_fit
    )
  })

  oriented_loading_df <- reactive({
    orient_loading_data(
      analysis()$loading_df,
      pc_sign_vector(input),
      pca_fit = analysis()$pca$pca_fit
    )
  })

  for (i in seq_len(6)) {
    local({
      pc_index <- i
      pc_id <- paste0("PC", pc_index)
      output[[paste0("pc", pc_index, "_spatial")]] <- renderPlot({
        req(pc_index <= length(analysis()$pca$var_explained))
        req(pc_id %in% names(oriented_spatial_df()))
        plot_pc_spatial(
          oriented_spatial_df(),
          pc_id,
          analysis()$pca$var_explained[[pc_index]],
          pc_id,
          analysis()$spatial_inputs$category_polygon_sf
        ) +
          labs(title = NULL) +
          theme(
            plot.title = element_blank(),
            plot.margin = margin(4, 4, 4, 4)
          )
      })
    })
  }

  output$pc_spatial_outputs <- renderUI({
    pc_ids <- selected_pc_ids(input, analysis())
    if (length(pc_ids) == 0) {
      return(app_panel_card(
        "No PCs Selected",
        class = "setup-empty",
        tags$p("Select one or more PCs in the sidebar to display spatial score maps.")
      ))
    }
    n_plots <- length(pc_ids)
    plot_ids <- paste0(tolower(pc_ids), "_spatial")

    tags$div(
      class = paste("pca-map-grid", paste0("pc-count-", n_plots)),
      lapply(seq_len(n_plots), function(i) {
        pc_index <- as.integer(sub("^PC", "", pc_ids[[i]]))
        orientation <- if (identical(input[[paste0("pc", pc_index, "_orientation")]], "reverse")) "Reversed" else "Keep"
        app_panel_card(
          pc_ids[[i]],
          class = "pca-map-card",
          tags$p(
            class = "pca-card-meta",
            paste0(sprintf("%.2f", 100 * analysis()$pca$var_explained[[pc_index]]), "% variance explained · Orientation: ", orientation)
          ),
          plotOutput(plot_ids[[i]], height = 330)
        )
      })
    )
  })

  output$pca_detail_meta <- renderUI({
    pc <- input$loading_pc %||% "PC1"
    pc_index <- as.integer(sub("^PC", "", pc))
    req(!is.na(pc_index), pc_index <= length(analysis()$pca$var_explained))
    orientation <- if (identical(input[[paste0("pc", pc_index, "_orientation")]], "reverse")) "Reversed" else "Keep"
    tags$p(
      class = "pca-card-meta",
      paste0(sprintf("%.2f", 100 * analysis()$pca$var_explained[[pc_index]]), "% variance explained · Orientation: ", orientation)
    )
  })

  output$pca_detail_spatial <- renderPlot({
    pc <- input$loading_pc %||% "PC1"
    pc_index <- as.integer(sub("^PC", "", pc))
    req(!is.na(pc_index), pc_index <= length(analysis()$pca$var_explained))
    req(pc %in% names(oriented_spatial_df()))

    plot_pc_spatial(
      oriented_spatial_df(),
      pc,
      analysis()$pca$var_explained[[pc_index]],
      pc,
      analysis()$spatial_inputs$category_polygon_sf
    ) +
      labs(title = NULL) +
      theme(
        plot.title = element_blank(),
        plot.margin = margin(4, 4, 4, 4)
      )
  })

  output$loading_pc_control <- renderUI({
    candidates <- pc_candidate_ids(analysis(), max_pcs = 6)
    if (length(candidates) == 0) {
      return(tags$p(class = "setup-meta", "No PCs are available to inspect."))
    }
    current <- input$loading_pc %||% candidates[[1]]
    if (!current %in% candidates) {
      current <- candidates[[1]]
    }
    selectInput("loading_pc", "PC to inspect", choices = candidates, selected = current)
  })

  loading_bar_panel <- function(title, loading_df, direction) {
    values <- abs(loading_df$loading)
    max_value <- max(values, na.rm = TRUE)
    if (!is.finite(max_value) || max_value == 0) {
      max_value <- 1
    }

	    tags$section(
	      class = paste("loading-panel", direction),
      h5(title),
      lapply(seq_len(nrow(loading_df)), function(i) {
        width <- 100 * abs(loading_df$loading[[i]]) / max_value
        tags$div(
          class = "loading-bar-row",
          tags$div(
            class = "loading-bar-head",
            tags$span(class = "loading-gene", loading_df$gene[[i]]),
            tags$span(class = "loading-value", formatC(loading_df$loading[[i]], format = "f", digits = 5))
          ),
          tags$div(
            class = "loading-bar-track",
            tags$div(
              class = paste("loading-bar-fill", direction),
              style = paste0("width:", sprintf("%.1f", width), "%;")
            )
          )
        )
      })
    )
  }

  output$loading_results <- renderUI({
    req(input$loading_pc, input$loading_n, input$loading_sign)
    n <- as.integer(input$loading_n)
    show_sign <- input$loading_sign %||% "both"
    pos <- top_loading_table(oriented_loading_df(), input$loading_pc, "positive", n = n)
    neg <- top_loading_table(oriented_loading_df(), input$loading_pc, "negative", n = n)

    if (identical(show_sign, "positive")) {
      return(tags$div(
        class = "loading-result-grid single",
        loading_bar_panel("Top Positive Loadings", pos, "positive")
      ))
    }
    if (identical(show_sign, "negative")) {
      return(tags$div(
        class = "loading-result-grid single",
        loading_bar_panel("Top Negative Loadings", neg, "negative")
      ))
    }

    tags$div(
      class = "loading-result-grid",
      loading_bar_panel("Top Positive Loadings", pos, "positive"),
      loading_bar_panel("Top Negative Loadings", neg, "negative")
    )
  })

  mediation_summary <- reactive({
    result <- frequentist_current_result()
    validate(need(!is.null(result), frequentist_unavailable_message()))
    result$result$fit$summary |>
      transmute(
        contrast,
        TE,
        direct_effect = NDE,
        total_indirect_effect = NIE,
        PM
      ) |>
      filter_by_mediation_inputs(
        contrast = input$mediation_contrast_filter %||% "all",
        mediator = "all"
      )
	  })

	  mediation_path <- reactive({
	    result <- frequentist_current_result()
	    validate(need(!is.null(result), frequentist_unavailable_message()))
	    result$result$fit$path |>
	      transmute(
        contrast,
        mediator,
        alpha = alpha_X_to_M,
        beta = beta_M_to_Y,
        indirect_component
	      ) |>
	      mutate(
	        mediator_id = mediator,
	        mediator = mediator_label(mediator)
	      ) |>
	      filter_by_mediation_inputs(
	        contrast = input$mediation_contrast_filter %||% "all",
	        mediator = input$mediation_mediator_filter %||% "all"
	      )
	  })

	  bayesian_mediation_summary <- reactive({
	    artifact <- bayesian_current_artifact()
	    validate(need(!is.null(artifact), bayesian_unavailable_message()))
	    get_bayesian_mediation_summary(artifact) |>
	      filter_by_mediation_inputs(
	        contrast = input$mediation_contrast_filter %||% "all",
	        mediator = "all"
	      )
	  })

	  bayesian_mediation_path <- reactive({
	    artifact <- bayesian_current_artifact()
	    validate(need(!is.null(artifact), bayesian_unavailable_message()))
	    path <- get_bayesian_mediation_path(artifact)
	    if (isTRUE(bayesian_artifact_available())) {
	      path <- orient_bayesian_path(path, mediation_sign_vector_confirmed())
	    } else {
	      path <- path |>
	        mutate(
	          mediator_id = mediator,
	          mediator = mediator_label(mediator)
	        )
	    }
	    path |>
	      filter_by_mediation_inputs(
	        contrast = input$mediation_contrast_filter %||% "all",
	        mediator = input$mediation_mediator_filter %||% "all"
	      )
	  })

	  bootstrap_cache_key <- function(freq, B) {
	    paste(freq$key, paste0("B=", as.integer(B)), sep = ":")
	  }

	  compute_bootstrap_result <- function(freq, B) {
	    boot <- withProgress(message = "Running bootstrap confidence intervals", value = 0, {
	      bootstrap_frequentist_dynamic(
	        frequentist_result = freq$result,
	        B = B,
	        seed = 123,
	        progress_callback = function(done, total) {
	          incProgress(
	            amount = 1 / total,
	            detail = paste0(done, " / ", total, " resamples")
	          )
	        }
	      )
	    })
	    list(
	      key = freq$key,
	      B = B,
	      ci_method = "empirical percentile quantiles at 2.5% and 97.5%",
	      ci = boot$ci,
	      bootstrap = boot$bootstrap
	    )
	  }

	  bootstrap_result <- reactive({
	    freq <- frequentist_current_result()
	    if (is.null(freq)) {
	      return(NULL)
	    }
	    B <- as.integer(input$bootstrap_B %||% 1000)
	    validate(need(B > 0, "Bootstrap replicates must be positive."))
	    key <- bootstrap_cache_key(freq, B)
	    cache <- frequentist_bootstrap_cache()
	    if (!is.null(cache[[key]])) {
	      return(cache[[key]])
	    }
	    result <- compute_bootstrap_result(freq, B)
	    cache[[key]] <- result
	    frequentist_bootstrap_cache(cache)
	    result
	  })

	  bootstrap_result_or_null <- reactive({
	    tryCatch(
	      bootstrap_result(),
	      shiny.silent.error = function(e) NULL,
      validation = function(e) NULL,
      error = function(e) NULL
    )
  })

  current_bootstrap_summary_ci <- reactive({
    boot <- bootstrap_result_or_null()
    if (is.null(boot) || !identical(boot$key, frequentist_current_cache_key())) {
      return(NULL)
    }
    boot$ci$summary_ci
  })

  current_bootstrap_path_ci <- reactive({
    boot <- bootstrap_result_or_null()
    if (is.null(boot) || !identical(boot$key, frequentist_current_cache_key())) {
      return(NULL)
    }
    boot$ci$path_ci
  })

  mediation_path_with_bootstrap_ci <- reactive({
    rows <- mediation_path()
    ci <- current_bootstrap_path_ci()
    if (is.null(ci)) {
      return(rows)
    }
    rows |>
      left_join(
        ci |>
          select(
            contrast,
            mediator,
            alpha_lwr,
            alpha_upr,
            beta_lwr,
            beta_upr,
            indirect_lwr,
            indirect_upr
          ),
        by = c("contrast", "mediator")
      )
  })

  mediation_interval_text <- function(lwr, upr, percent = FALSE) {
    if (is.null(lwr) || is.null(upr) || is.na(lwr) || is.na(upr)) {
      return(NULL)
    }
    if (isTRUE(percent)) {
      return(paste0("[", format_number(100 * lwr, 1), "%, ", format_number(100 * upr, 1), "%]"))
    }
    paste0("[", format_number(lwr, 3), ", ", format_number(upr, 3), "]")
  }

  mediation_value_text <- function(x, percent = FALSE) {
    if (isTRUE(percent)) {
      return(paste0(format_number(100 * x, 1), "%"))
    }
    format_number(x, 3)
  }

  mediation_probability_text <- function(p) {
    if (is.null(p) || is.na(p)) {
      return(NULL)
    }
    paste0("Pr(effect > 0): ", format_number(100 * p, 1), "%")
  }

  mediation_effect_card <- function(label, value, lwr = NA_real_, upr = NA_real_, probability = NULL, percent = FALSE) {
    interval <- mediation_interval_text(lwr, upr, percent = percent)
    tags$div(
      class = "effect-card",
      tags$div(class = "effect-label", label),
      tags$div(class = "effect-value", mediation_value_text(value, percent = percent)),
      if (!is.null(interval)) tags$div(class = "effect-interval", interval),
      if (!is.null(probability)) tags$div(class = "effect-evidence", mediation_probability_text(probability))
    )
  }

  mediation_summary_panel <- function(contrast, cards) {
    app_panel_card(
      contrast,
      class = "mediation-contrast-panel",
      tags$div(class = "effect-card-grid", cards)
    )
  }

  output$mediation_summary_cards <- renderUI({
    framework <- input$analysis_framework %||% "bayesian"
    if (identical(framework, "bayesian")) {
      rows <- bayesian_mediation_summary()
      return(tagList(lapply(seq_len(nrow(rows)), function(i) {
        row <- rows[i, ]
        mediation_summary_panel(
          row$contrast,
          list(
            mediation_effect_card("Total Effect", row$TE_est, row$TE_lwr, row$TE_upr, row$TE_Pr_gt_0),
            mediation_effect_card("Direct Effect", row$direct_est, row$direct_lwr, row$direct_upr, row$direct_Pr_gt_0),
            mediation_effect_card("Indirect Effect", row$total_IIE_est, row$total_IIE_lwr, row$total_IIE_upr, row$total_IIE_Pr_gt_0),
            mediation_effect_card("Proportion Mediated", row$PM_est, row$PM_lwr, row$PM_upr, percent = TRUE)
          )
        )
      })))
    }

    rows <- mediation_summary()
    ci <- current_bootstrap_summary_ci()
    if (!is.null(ci)) {
      rows <- rows |>
        left_join(
          ci |> select(contrast, TE_lwr, TE_upr, NDE_lwr, NDE_upr, NIE_lwr, NIE_upr, PM_lwr, PM_upr),
          by = "contrast"
        )
    } else {
      rows <- rows |>
        mutate(
          TE_lwr = NA_real_, TE_upr = NA_real_,
          NDE_lwr = NA_real_, NDE_upr = NA_real_,
          NIE_lwr = NA_real_, NIE_upr = NA_real_,
          PM_lwr = NA_real_, PM_upr = NA_real_
        )
    }

    tagList(lapply(seq_len(nrow(rows)), function(i) {
      row <- rows[i, ]
      mediation_summary_panel(
        row$contrast,
        list(
          mediation_effect_card("Total Effect", row$TE, row$TE_lwr, row$TE_upr),
          mediation_effect_card("Direct Effect", row$direct_effect, row$NDE_lwr, row$NDE_upr),
          mediation_effect_card("Indirect Effect", row$total_indirect_effect, row$NIE_lwr, row$NIE_upr),
          mediation_effect_card("Proportion Mediated", row$PM, row$PM_lwr, row$PM_upr, percent = TRUE)
        )
      )
    }))
  })

  bayesian_path_evidence <- function(pr_gt, pr_lt) {
    if (is.na(pr_gt) && is.na(pr_lt)) {
      return("")
    }
    if (!is.na(pr_gt) && (is.na(pr_lt) || pr_gt >= pr_lt)) {
      return(paste0("Pr(positive) = ", format_number(100 * pr_gt, 1), "%"))
    }
    paste0("Pr(negative) = ", format_number(100 * pr_lt, 1), "%")
  }

  frequentist_path_cell <- function(est, lwr = NA_real_, upr = NA_real_, digits = 4) {
    if (is.na(lwr) || is.na(upr)) {
      return(format_number(est, digits = digits))
    }
    format_estimate_interval(est, lwr, upr, digits = digits)
  }

  mediation_path_table_ui <- function(df, framework) {
    if (nrow(df) == 0) {
      return(tags$p(class = "setup-meta", "No path rows match the current filters."))
    }
	    header <- if (identical(framework, "bayesian")) {
	      c("Contrast", "Mediator", "Exposure -> Mediator", "Mediator -> Outcome", "Indirect Effect", "Posterior Direction")
	    } else {
	      c("Contrast", "Mediator", "Exposure -> Mediator", "Mediator -> Outcome", "Indirect Effect")
	    }
    rows <- lapply(seq_len(nrow(df)), function(i) {
      if (identical(framework, "bayesian")) {
        tags$tr(
          tags$td(df$contrast[[i]]),
          tags$td(df$mediator[[i]]),
          tags$td(format_estimate_interval(df$alpha_est[[i]], df$alpha_lwr[[i]], df$alpha_upr[[i]])),
          tags$td(format_estimate_interval(df$beta_est[[i]], df$beta_lwr[[i]], df$beta_upr[[i]])),
          tags$td(format_estimate_interval(df$IIE_est[[i]], df$IIE_lwr[[i]], df$IIE_upr[[i]])),
          tags$td(bayesian_path_evidence(df$IIE_Pr_gt_0[[i]], df$IIE_Pr_lt_0[[i]]))
        )
      } else {
        value_or_na <- function(column) {
          if (!column %in% names(df)) {
            return(NA_real_)
          }
          df[[column]][[i]]
        }
        tags$tr(
          tags$td(df$contrast[[i]]),
          tags$td(df$mediator[[i]]),
          tags$td(frequentist_path_cell(df$alpha[[i]], value_or_na("alpha_lwr"), value_or_na("alpha_upr"), digits = 4)),
          tags$td(frequentist_path_cell(df$beta[[i]], value_or_na("beta_lwr"), value_or_na("beta_upr"), digits = 4)),
          tags$td(frequentist_path_cell(df$indirect_component[[i]], value_or_na("indirect_lwr"), value_or_na("indirect_upr"), digits = 4))
        )
      }
    })
    tags$table(
      class = "mediation-path-table",
      tags$thead(tags$tr(lapply(header, tags$th))),
      tags$tbody(rows)
    )
  }

  output$mediation_path_details <- renderUI({
    framework <- input$analysis_framework %||% "bayesian"
    if (identical(framework, "bayesian")) {
      return(mediation_path_table_ui(bayesian_mediation_path(), "bayesian"))
    }
    mediation_path_table_ui(mediation_path_with_bootstrap_ci(), "frequentist")
  })

	  output$mediation_main_content <- renderUI({
	    framework <- input$analysis_framework %||% "bayesian"
	    result_available <- if (identical(framework, "bayesian")) {
	      !is.null(bayesian_current_artifact())
	    } else {
	      !is.null(frequentist_current_result())
	    }
	    if (!isTRUE(result_available)) {
	      button_id <- if (identical(framework, "bayesian")) "run_bayesian_mediation" else "run_frequentist_mediation"
	      button_label <- if (identical(framework, "bayesian")) "Run Bayesian Mediation" else "Run Frequentist Mediation"
	      running <- if (identical(framework, "bayesian")) isTRUE(bayesian_fit_running()) else isTRUE(frequentist_fit_running())
	      hint <- if (identical(framework, "bayesian")) bayesian_unavailable_message() else frequentist_unavailable_message()
		      return(app_panel_card(
		          "No mediation results yet",
		          class = "setup-empty",
		          tags$p(paste0(
		            "Run ",
		            if (identical(framework, "bayesian")) "Bayesian" else "Frequentist",
		            " Mediation using the controls on the left to estimate mediation effects for the current data and selected PCs."
		          )),
		          tags$p(class = "setup-meta", hint),
		          tags$div(
		            class = "mediation-empty-action",
		            actionButton(button_id, if (running) paste(button_label, "running") else button_label, class = "btn-primary")
		          )
		      ))
		    }

    framework_label <- if (identical(framework, "bayesian")) "Bayesian posterior estimates" else "Frequentist estimates"
    interval_label <- if (identical(framework, "bayesian")) {
      "Posterior mean with 95% credible interval."
    } else if (is.null(current_bootstrap_summary_ci())) {
	      "Point estimates shown while bootstrap confidence intervals are prepared."
    } else {
      "Estimate with 95% bootstrap confidence interval."
    }
    contrast <- input$mediation_contrast_filter %||% "all"
    contrast_context <- if (identical(contrast, "all")) "All contrasts" else contrast
    plot_id <- if (identical(framework, "bayesian")) "bayesian_indirect_plot" else "indirect_plot"

	    tagList(
	      tags$div(
	        class = "mediation-header",
	        tags$div(
	          h3("Mediation Results"),
	          tags$p(class = "setup-meta", framework_label)
	        ),
	        tags$div(
	          class = "mediation-header-actions",
	          tags$div(class = "mediation-context", contrast_context),
	          actionButton("continue_to_sensitivity", "Sensitivity", class = "btn-default")
	        )
	      ),
      tags$div(
        class = "mediation-result-grid",
        div(
          class = "mediation-view-tabs",
          tabsetPanel(
            id = "mediation_result_view",
            type = "tabs",
            tabPanel(
              "Overall Effects",
              uiOutput("mediation_summary_cards")
            ),
            tabPanel(
              "Indirect Effects",
              app_panel_card(
                "Indirect Effects",
                tags$p(class = "setup-meta", interval_label),
                plotOutput(plot_id, height = 430)
              )
            ),
            tabPanel(
              "Path Details",
              app_panel_card(
                "Path Details",
                tags$p(
                  class = "setup-meta",
                  if (identical(framework, "bayesian")) {
                    "Path estimates are posterior means with 95% credible intervals."
                  } else if (is.null(current_bootstrap_path_ci())) {
                    "Path estimates use the current fitted mediation model."
                  } else {
                    "Path estimates with 95% bootstrap confidence intervals."
                  }
                ),
                uiOutput("mediation_path_details")
              )
            )
          )
        )
      )
    )
  })

  order_mediators_for_forest_plot <- function(df) {
    mediators_in_order <- df$mediator |>
      unique() |>
      sort()
    df$mediator <- factor(df$mediator, levels = rev(mediators_in_order))
    df
  }

	  output$bayesian_mediation_summary <- renderDT({
	    bayesian_mediation_summary() |>
	      present_bayesian_mediation_summary() |>
	      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "tip", scrollX = FALSE))
	  }, server = FALSE)

	  output$bayesian_mediation_probabilities <- renderDT({
	    bayesian_mediation_summary() |>
	      present_bayesian_mediation_probabilities() |>
	      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "t", scrollX = FALSE))
	  }, server = FALSE)

	  output$bayesian_mediation_path <- renderDT({
	    bayesian_mediation_path() |>
	      present_bayesian_mediation_path() |>
	      datatable(rownames = FALSE, options = list(pageLength = 10, dom = "tip", scrollX = FALSE))
	  }, server = FALSE)

	  output$bayesian_indirect_plot <- renderPlot({
	    df <- bayesian_mediation_path()
      validate(need(nrow(df) > 0, "No mediator-specific effects match the current filters."))
      df <- order_mediators_for_forest_plot(df)
      dodge <- position_dodge(width = 0.55)
	    ggplot(df, aes(x = IIE_est, y = mediator, color = contrast)) +
	      geom_vline(xintercept = 0, linewidth = 0.35, color = "grey60") +
	      geom_errorbarh(
	        aes(xmin = IIE_lwr, xmax = IIE_upr),
	        position = dodge,
	        height = 0.18,
	        linewidth = 0.55
	      ) +
	      geom_point(position = dodge, size = 2.3) +
	      labs(x = "Indirect Effect, posterior mean with 95% CrI", y = "Mediator", color = "Contrast") +
	      theme_bw(base_size = 12) +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = if (dplyr::n_distinct(df$contrast) > 1) "bottom" else "none"
        )
	  })

	  output$mediation_summary <- renderDT({
	    mediation_summary() |>
      present_frequentist_mediation_summary() |>
      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "tip", scrollX = FALSE))
  }, server = FALSE)

  output$mediation_path <- renderDT({
    mediation_path() |>
      present_frequentist_mediation_path() |>
      datatable(rownames = FALSE, options = list(pageLength = 10, dom = "tip", scrollX = FALSE))
  }, server = FALSE)

  output$indirect_plot <- renderPlot({
    df <- mediation_path_with_bootstrap_ci()
    validate(need(nrow(df) > 0, "No mediator-specific effects match the current filters."))
    df <- order_mediators_for_forest_plot(df)
    dodge <- position_dodge(width = 0.55)
    p <- ggplot(df, aes(x = indirect_component, y = mediator, color = contrast)) +
      geom_vline(xintercept = 0, linewidth = 0.35, color = "grey60") +
      labs(
        x = if (!is.null(current_bootstrap_path_ci())) "Indirect Effect, estimate with 95% bootstrap CI" else "Indirect Effect",
        y = "Mediator",
        color = "Contrast"
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = if (dplyr::n_distinct(df$contrast) > 1) "bottom" else "none"
      )
    if (all(c("indirect_lwr", "indirect_upr") %in% names(df)) && any(!is.na(df$indirect_lwr) & !is.na(df$indirect_upr))) {
      p <- p +
        geom_errorbarh(
          aes(xmin = indirect_lwr, xmax = indirect_upr),
          position = dodge,
          height = 0.18,
          linewidth = 0.55
        )
    }
    p + geom_point(position = dodge, size = 2.4)
  })

	  output$bootstrap_status <- renderUI({
	    boot <- bootstrap_result_or_null()
	    if (is.null(boot)) {
	      return(tags$div(class = "setup-meta", "Shown automatically after frequentist results are available."))
	    } else if (!identical(boot$key, frequentist_current_cache_key())) {
	      return(tags$div(class = "setup-meta", style = "color:#8a5a00;", "Updating for the current setup."))
	    } else {
	      tags$div(class = "setup-meta", style = "color:#1b6e3c;", paste("Automatic 95% CI with B =", boot$B))
	    }
	  })

  output$bootstrap_summary_ci <- renderDT({
    req(bootstrap_result())
    validate(need(identical(bootstrap_result()$key, frequentist_current_cache_key()), "Bootstrap CI is updating for the current setup."))
    bootstrap_result()$ci$summary_ci |>
      dplyr::filter((input$mediation_contrast_filter %||% "all") == "all" | contrast == (input$mediation_contrast_filter %||% "all")) |>
      present_bootstrap_summary_ci(contrast_labels = NULL) |>
      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "tip", scrollX = FALSE))
  }, server = FALSE)

	  output$bootstrap_path_ci <- renderDT({
	    req(bootstrap_result())
	    validate(need(identical(bootstrap_result()$key, frequentist_current_cache_key()), "Bootstrap CI is updating for the current setup."))
	    bootstrap_result()$ci$path_ci |>
	      dplyr::filter(
	        (input$mediation_contrast_filter %||% "all") == "all" | contrast == (input$mediation_contrast_filter %||% "all"),
	        (input$mediation_mediator_filter %||% "all") == "all" | mediator == (input$mediation_mediator_filter %||% "all")
	      ) |>
	      present_bootstrap_path_ci(contrast_labels = NULL) |>
	      datatable(rownames = FALSE, options = list(pageLength = 10, dom = "tip", scrollX = FALSE))
	  }, server = FALSE)

	  output$sensitivity_main_content <- renderUI({
	    framework <- input$sensitivity_framework %||% "bayesian"
	    scenario <- input$sensitivity_scenario %||% "common"
	    rho <- input$sensitivity_active_rho %||% 0
	    result_available <- if (identical(framework, "bayesian")) {
	      !is.null(bayesian_current_artifact())
	    } else {
	      !is.null(frequentist_current_result())
	    }

	    if (!isTRUE(result_available)) {
	      return(app_panel_card(
	          if (identical(framework, "bayesian")) "Bayesian mediation is not fitted yet" else "Frequentist mediation is not fitted yet",
	          class = "setup-empty",
	          tags$p(
	            if (identical(framework, "bayesian")) {
	              "Run Bayesian Mediation first to enable Bayesian sensitivity analysis for the current data and selected PCs."
	            } else {
	              "Run Frequentist Mediation first to enable Frequentist sensitivity analysis for the current data and selected PCs."
	            }
	          ),
	          actionButton("go_to_mediation_from_sensitivity", "Go to Mediation", class = "btn-primary")
	      ))
	    }

    framework_label <- if (identical(framework, "bayesian")) {
      "Mediator-outcome confounding sensitivity with posterior mean and 95% credible interval."
    } else {
      "Frequentist point-estimate sensitivity analysis."
    }
    scenario_label <- if (identical(scenario, "common")) {
      "Common correlation"
    } else {
      paste("One mediator:", mediator_label(input$sensitivity_one_mediator %||% ""))
    }
    primary_plot_id <- if (identical(scenario, "common")) {
      if (identical(framework, "bayesian")) "bayesian_sensitivity_curve" else "sensitivity_curve"
    } else {
      if (identical(framework, "bayesian")) "bayesian_sensitivity_target_curve" else "sensitivity_target_curve"
    }
    primary_plot_title <- if (identical(scenario, "common")) {
      "Total Indirect Effect"
    } else {
      paste(mediator_label(input$sensitivity_one_mediator %||% ""), "Indirect Effect")
    }
	    path_plot_id <- if (identical(scenario, "common")) {
	      if (identical(framework, "bayesian")) "bayesian_sensitivity_path_curve" else "sensitivity_path_curve"
	    } else {
	      NULL
	    }
    warning_id <- if (identical(framework, "bayesian")) "bayesian_sensitivity_warnings" else "frequentist_sensitivity_warnings"

	    tagList(
	      tags$div(
	        class = "sensitivity-header",
        tags$div(
          h3("Sensitivity Analysis"),
          tags$p(class = "setup-meta", framework_label)
        ),
        tags$div(class = "sensitivity-context", paste0(scenario_label, " · rho = ", format_number(rho, digits = 2)))
      ),
      tags$div(
        class = "sensitivity-result-grid",
        div(
          class = "sensitivity-view-tabs",
          tabsetPanel(
            id = "sensitivity_result_view",
            type = "tabs",
            tabPanel(
              primary_plot_title,
              app_panel_card(
                "Zero Crossing",
                uiOutput("sensitivity_tipping_summary")
              ),
              app_panel_card(
                primary_plot_title,
                plotOutput(primary_plot_id, height = 430)
              )
            ),
            tabPanel(
              "Selected rho Decomposition",
              app_panel_card(
                "Selected rho Decomposition",
                tags$p(class = "setup-meta", if (identical(framework, "bayesian")) "Posterior mean with 95% credible interval." else "Point estimates at the selected rho."),
                uiOutput("sensitivity_selected_decomposition")
              ),
              uiOutput(warning_id)
            ),
            if (identical(scenario, "common")) {
              tabPanel(
                "Mediator Components",
                app_panel_card(
                  "Mediator Components",
                  tags$p(
                    class = "setup-meta",
                    "The same residual correlation is applied to all mediator components."
                  ),
                  plotOutput(path_plot_id, height = 430)
                )
              )
            },
	            if (identical(scenario, "one_at_a_time")) {
	              tabPanel(
	                "All Mediators Comparison",
	                app_panel_card(
	                  "All Mediators Comparison",
	                  tags$p(
	                    class = "setup-meta",
	                    "Each curve varies one mediator at a time and shows that mediator's own indirect effect."
	                  ),
	                  plotOutput(
	                    if (identical(framework, "bayesian")) "bayesian_sensitivity_all_mediator_comparison" else "sensitivity_all_mediator_comparison",
	                    height = 430
	                  )
	                )
	              )
	            }
	          )
	        )
	      )
    )
  })

	  bayesian_sensitivity_selected <- reactive({
	    req(identical(input$sensitivity_framework, "bayesian"))
	    artifact <- bayesian_current_artifact()
	    validate(need(!is.null(artifact), "Run Bayesian Mediation first for the current data and selected PCs."))
      result <- bayesian_current_result()
	    scenario <- input$sensitivity_scenario %||% "common"
	    rho <- input$sensitivity_active_rho %||% 0

	    if (identical(scenario, "common")) {
	      bayesian_sensitivity_selected_cache(
	        artifact = artifact,
          artifact_key = result$key,
	        rho = rho,
	        mode = "common"
	      )
	    } else {
	      mediator <- input$sensitivity_one_mediator %||% bayes_A_bundle_mediators(artifact$sensitivity_A$residual_draws)[[1]]
	      bayesian_sensitivity_selected_cache(
	        artifact = artifact,
          artifact_key = result$key,
	        rho = rho,
	        mode = "one_at_a_time",
	        mediator = mediator
	      )
	    }
	  })

  sensitivity_path_selected <- reactive({
    req(identical(input$sensitivity_framework, "frequentist"))
    req(input$sensitivity_scenario, input$sensitivity_active_rho)

    if (identical(input$sensitivity_scenario, "common")) {
      sensitivity()$common$sensitivity_path |>
        filter(valid_parameter, abs(rho - input$sensitivity_active_rho) < 1e-10) |>
        mutate(contrast = contrast_label_current(exposure), mediator = mediator_label(mediator)) |>
        filter_sensitivity_contrast() |>
        transmute(
          contrast,
          rho,
          mediator,
          indirect_est = indirect_adjusted
        )
    } else {
      req(input$sensitivity_one_mediator)
      sensitivity()$one_at_a_time$sensitivity_path |>
        filter(
          valid_parameter,
          varied_mediator == input$sensitivity_one_mediator,
          mediator == input$sensitivity_one_mediator,
          abs(rho - input$sensitivity_active_rho) < 1e-10
        ) |>
        mutate(contrast = contrast_label_current(exposure), mediator = mediator_label(mediator)) |>
        filter_sensitivity_contrast() |>
        transmute(
          contrast,
          rho,
          mediator,
          indirect_est = indirect_adjusted
        )
    }
  })

  sensitivity_value_cell <- function(est, lwr = NA_real_, upr = NA_real_, percent = FALSE) {
    if (isTRUE(percent)) {
      est_txt <- paste0(format_number(100 * est, 1), "%")
      if (is.na(lwr) || is.na(upr)) {
        return(est_txt)
      }
      return(paste0(est_txt, " [", format_number(100 * lwr, 1), "%, ", format_number(100 * upr, 1), "%]"))
    }
    format_estimate_interval(est, lwr, upr)
  }

  sensitivity_mediator_table <- function(path_df, framework = c("bayesian", "frequentist")) {
    framework <- match.arg(framework)
    if (nrow(path_df) == 0) {
      return(tags$p(class = "setup-meta", "No mediator-specific values are feasible at the selected rho."))
    }
    rows <- lapply(seq_len(nrow(path_df)), function(i) {
      if (identical(framework, "bayesian")) {
        tags$tr(
          tags$td(mediator_label(path_df$mediator[[i]])),
          tags$td(format_estimate_interval(path_df$IIE_est[[i]], path_df$IIE_lwr[[i]], path_df$IIE_upr[[i]]))
        )
      } else {
        tags$tr(
          tags$td(path_df$mediator[[i]]),
          tags$td(format_number(path_df$indirect_est[[i]]))
        )
      }
    })
    tags$table(
      class = "sensitivity-mini-table",
      tags$thead(tags$tr(tags$th("Mediator"), tags$th("Indirect Effect"))),
      tags$tbody(rows)
    )
  }

  first_or_empty <- function(x) {
    if (length(x) == 0 || is.null(x[[1]]) || is.na(x[[1]])) "" else x[[1]]
  }

  zero_crossing_rho <- function(rho, value) {
    ok <- is.finite(rho) & is.finite(value)
    rho <- rho[ok]
    value <- value[ok]
    if (length(rho) == 0) {
      return(NA_real_)
    }
    ord <- order(rho)
    rho <- rho[ord]
    value <- value[ord]

    exact <- which(abs(value) < .Machine$double.eps^0.5)
    if (length(exact) > 0) {
      return(rho[[exact[[1]]]])
    }
    if (length(rho) < 2) {
      return(NA_real_)
    }

    changes <- which(value[-length(value)] * value[-1] < 0)
    if (length(changes) == 0) {
      return(NA_real_)
    }
    i <- changes[[1]]
    rho[[i]] - value[[i]] * (rho[[i + 1]] - rho[[i]]) / (value[[i + 1]] - value[[i]])
  }

  sensitivity_primary_curve_data <- reactive({
    framework <- input$sensitivity_framework %||% "bayesian"
    scenario <- input$sensitivity_scenario %||% "common"

    if (identical(framework, "bayesian")) {
      if (identical(scenario, "common")) {
        return(bayesian_sensitivity_curve()$summary |>
          filter_sensitivity_contrast() |>
          transmute(contrast, rho, value = total_IIE_est))
      }
      target <- input$sensitivity_one_mediator %||% ""
      return(bayesian_sensitivity_curve()$path |>
        dplyr::filter(mediator == target) |>
        filter_sensitivity_contrast() |>
        transmute(contrast, rho, value = IIE_est))
    }

    if (identical(scenario, "common")) {
      return(sensitivity()$common$sensitivity_summary |>
        filter(valid_parameter) |>
        mutate(contrast = contrast_label_current(exposure)) |>
        filter_sensitivity_contrast() |>
        transmute(contrast, rho, value = NIE_adjusted))
    }

    target <- input$sensitivity_one_mediator %||% ""
    sensitivity()$one_at_a_time$sensitivity_path |>
      filter(valid_parameter, varied_mediator == target, mediator == target) |>
      mutate(contrast = contrast_label_current(exposure)) |>
      filter_sensitivity_contrast() |>
      transmute(contrast, rho, value = indirect_adjusted)
  })

  output$sensitivity_tipping_summary <- renderUI({
    curve <- sensitivity_primary_curve_data()
    validate(need(nrow(curve) > 0, "No sensitivity curve rows match the current filters."))

    tipping <- curve |>
      group_by(contrast) |>
      summarise(
        rho_zero_crossing = zero_crossing_rho(rho, value),
        crossing_in_grid = is.finite(rho_zero_crossing),
        .groups = "drop"
      )

    scenario <- input$sensitivity_scenario %||% "common"
    quantity <- if (identical(scenario, "common")) {
      "NIE = 0"
    } else {
      paste(mediator_label(input$sensitivity_one_mediator %||% ""), "indirect effect = 0")
    }

    tags$div(
      tags$p(class = "setup-meta", paste("Zero crossing for the primary curve:", quantity)),
      tags$div(
        class = "sensitivity-tipping-grid",
        lapply(seq_len(nrow(tipping)), function(i) {
          tags$div(
            class = "sensitivity-tipping-card",
            tags$div(class = "sensitivity-tipping-label", tipping$contrast[[i]]),
            tags$div(
              class = "sensitivity-tipping-value",
              if (isTRUE(tipping$crossing_in_grid[[i]])) {
                paste0("rho = ", format_number(tipping$rho_zero_crossing[[i]], digits = 3))
              } else {
                "No crossing in grid"
              }
            )
          )
        })
      )
    )
  })

  output$sensitivity_selected_decomposition <- renderUI({
    framework <- input$sensitivity_framework %||% "bayesian"
    scenario <- input$sensitivity_scenario %||% "common"

    if (identical(framework, "bayesian")) {
      selected <- bayesian_sensitivity_selected()
      summary <- filter_sensitivity_contrast(selected$summary)
      path <- filter_sensitivity_contrast(selected$path)
      if (identical(scenario, "one_at_a_time")) {
        target <- input$sensitivity_one_mediator %||% first_or_empty(selected$path$mediator_varied)
        path <- dplyr::filter(path, mediator == target)
      }

      return(tagList(lapply(seq_len(nrow(summary)), function(i) {
        row <- summary[i, ]
        path_rows <- dplyr::filter(path, contrast == row$contrast)
        app_panel_card(
          row$contrast,
          class = "sensitivity-decomp-panel",
          tags$div(
            class = "sensitivity-effect-grid",
            mediation_effect_card("Direct Effect", row$direct_est, row$direct_lwr, row$direct_upr),
            mediation_effect_card("Indirect Effect (Total)", row$total_IIE_est, row$total_IIE_lwr, row$total_IIE_upr),
            mediation_effect_card("Total Effect", row$TE_est, row$TE_lwr, row$TE_upr),
            mediation_effect_card("Proportion Mediated", row$PM_est, row$PM_lwr, row$PM_upr, percent = TRUE)
          ),
          h4(if (identical(scenario, "common")) "Mediator Components" else paste("Selected Mediator:", mediator_label(first_or_empty(path_rows$mediator)))),
          sensitivity_mediator_table(path_rows, framework = "bayesian")
        )
      })))
    }

    summary <- sensitivity_selected()
    path <- sensitivity_path_selected()
    tagList(lapply(seq_len(nrow(summary)), function(i) {
      row <- summary[i, ]
      path_rows <- dplyr::filter(path, contrast == row$contrast)
      app_panel_card(
        row$contrast,
        class = "sensitivity-decomp-panel",
        tags$div(
          class = "sensitivity-effect-grid",
          mediation_effect_card("Direct Effect", row$adjusted_direct_effect),
          mediation_effect_card("Indirect Effect (Total)", row$adjusted_total_indirect),
          mediation_effect_card("Total Effect", row$adjusted_total_effect),
          mediation_effect_card("Proportion Mediated", row$PM_adjusted, percent = TRUE)
          ),
        h4(if (identical(scenario, "common")) "Mediator Components" else paste("Selected Mediator:", first_or_empty(path_rows$mediator))),
        sensitivity_mediator_table(path_rows, framework = "frequentist")
      )
    }))
  })

	  bayesian_sensitivity_curve <- reactive({
	    req(identical(input$sensitivity_framework, "bayesian"))
	    artifact <- bayesian_current_artifact()
	    validate(need(!is.null(artifact), "Run Bayesian Mediation first for the current data and selected PCs."))
	    result <- bayesian_current_result()
	    scenario <- input$sensitivity_scenario %||% "common"

	    if (identical(scenario, "common")) {
	      bayesian_sensitivity_curve_cache(artifact = artifact, artifact_key = result$key, mode = "common")
	    } else {
	      mediator <- input$sensitivity_one_mediator %||% bayes_A_bundle_mediators(artifact$sensitivity_A$residual_draws)[[1]]
	      bayesian_sensitivity_curve_cache(artifact = artifact, artifact_key = result$key, mode = "one_at_a_time", mediator = mediator)
	    }
	  })

	  output$bayesian_sensitivity_summary <- renderDT({
	    bayesian_sensitivity_selected()$summary |>
	      present_bayesian_sensitivity_summary() |>
	      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "tip", scrollX = FALSE))
	  }, server = FALSE)

	  output$bayesian_sensitivity_path <- renderDT({
	    bayesian_sensitivity_selected()$path |>
	      present_bayesian_sensitivity_path() |>
	      datatable(rownames = FALSE, options = list(pageLength = 10, dom = "tip", scrollX = FALSE))
	  }, server = FALSE)

	  output$bayesian_sensitivity_warnings <- renderUI({
	    selected <- bayesian_sensitivity_selected()
	    warnings <- list()
	    feasible_fraction <- min(selected$summary$valid_draw_fraction, na.rm = TRUE)

	    if (is.finite(feasible_fraction) && feasible_fraction <= 0.25) {
	      warnings <- c(
	        warnings,
	        list(tags$p("This sensitivity setting is poorly feasible for the posterior covariance structure; do not interpret it as an ordinary valid result."))
	      )
	    } else if (is.finite(feasible_fraction) && feasible_fraction < 1) {
	      warnings <- c(
	        warnings,
	        list(tags$p("This sensitivity setting is not feasible for all posterior draws; interpret the result with caution."))
	      )
	    }

	    if (pm_warning_needed(selected$summary)) {
	      warnings <- c(
	        warnings,
	        list(tags$p("Proportion mediated is difficult to interpret at this sensitivity setting because the total effect is near zero or the direct and indirect effects have opposing signs."))
	      )
	    }

	    if (length(warnings) == 0) {
	      return(NULL)
	    }

	    tags$div(style = "color:#8a5a00;", warnings)
	  })

	  output$bayesian_sensitivity_curve <- renderPlot({
	    curve <- bayesian_sensitivity_curve()$summary |>
        filter_sensitivity_contrast()
	    active_rho <- input$sensitivity_active_rho %||% 0
      validate(need(nrow(curve) > 0, "No sensitivity curve rows match the current contrast filter."))

	    ggplot(curve, aes(x = rho, y = total_IIE_est, color = contrast, fill = contrast)) +
	      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
	      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
	      geom_vline(xintercept = active_rho, linetype = "dashed", color = "grey45") +
	      geom_ribbon(aes(ymin = total_IIE_lwr, ymax = total_IIE_upr), alpha = 0.12, color = NA) +
	      geom_line(linewidth = 0.9) +
	      geom_point(size = 1.4) +
	      labs(x = "Residual correlation (rho)", y = "Total Indirect Effect, posterior mean with 95% CrI", color = "Contrast", fill = "Contrast") +
	      theme_bw(base_size = 12)
	  })

    output$bayesian_sensitivity_target_curve <- renderPlot({
      req(identical(input$sensitivity_scenario, "one_at_a_time"))
      target <- input$sensitivity_one_mediator %||% ""
      curve <- bayesian_sensitivity_curve()$path |>
        dplyr::filter(mediator == target) |>
        filter_sensitivity_contrast()
      active_rho <- input$sensitivity_active_rho %||% 0
      validate(need(nrow(curve) > 0, "No target mediator sensitivity curve rows match the current filters."))

      ggplot(curve, aes(x = rho, y = IIE_est, color = contrast, fill = contrast)) +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
        geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
        geom_vline(xintercept = active_rho, linetype = "dashed", color = "grey45") +
        geom_ribbon(aes(ymin = IIE_lwr, ymax = IIE_upr), alpha = 0.12, color = NA) +
        geom_line(linewidth = 0.9) +
        geom_point(size = 1.4) +
        labs(
          x = paste("Residual correlation for", mediator_label(target)),
          y = paste(mediator_label(target), "Indirect Effect, posterior mean with 95% CrI"),
          color = "Contrast",
          fill = "Contrast"
        ) +
        theme_bw(base_size = 12)
    })

		  output$bayesian_sensitivity_path_curve <- renderPlot({
		    curve <- bayesian_sensitivity_curve()$path |>
	        filter_sensitivity_contrast()
	    selected_varied <- input$sensitivity_one_mediator %||% (unique(curve$mediator)[[1]] %||% "")
	    curve <- curve |>
	      mutate(
	        mediator_id = mediator,
	        mediator = mediator_label(mediator),
	        role = if_else(
	          identical(input$sensitivity_scenario, "one_at_a_time") & mediator_id == selected_varied,
	          "Varied mediator",
	          "Other mediators"
	        )
	      )
	    active_rho <- input$sensitivity_active_rho %||% 0
      validate(need(nrow(curve) > 0, "No mediator-specific sensitivity curve rows match the current filters."))

	    ggplot(curve, aes(x = rho, y = IIE_est, color = mediator, fill = mediator)) +
	      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
	      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
	      geom_vline(xintercept = active_rho, linetype = "dashed", color = "grey45") +
	      geom_ribbon(aes(ymin = IIE_lwr, ymax = IIE_upr), alpha = 0.10, color = NA) +
	      geom_line(aes(linetype = role, alpha = role), linewidth = 0.8) +
	      scale_alpha_manual(values = c("Varied mediator" = 1, "Other mediators" = 0.55), guide = "none") +
	      facet_wrap(~ contrast, scales = "free_y") +
		      labs(x = "Residual correlation (rho)", y = "Indirect Effect, posterior mean with 95% CrI", color = "Mediator", fill = "Mediator", linetype = "Role") +
		      theme_bw(base_size = 12)
		  })

	  output$bayesian_sensitivity_all_mediator_comparison <- renderPlot({
	    req(identical(input$sensitivity_framework, "bayesian"))
	    req(identical(input$sensitivity_scenario, "one_at_a_time"))
	    artifact <- bayesian_current_artifact()
	    validate(need(!is.null(artifact), "Run Bayesian Mediation first for the current data and selected PCs."))
	    result <- bayesian_current_result()
	    curve <- bayesian_sensitivity_all_mediator_curve_cache(
	      artifact = artifact,
	      artifact_key = result$key
	    )$path |>
	      dplyr::filter(mediator == mediator_varied) |>
	      filter_sensitivity_contrast() |>
	      mutate(mediator = mediator_label(mediator))
	    validate(need(nrow(curve) > 0, "No mediator comparison rows match the current contrast filter."))

	    ggplot(curve, aes(x = rho, y = IIE_est, color = mediator, fill = mediator)) +
	      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
	      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
	      geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45") +
	      geom_ribbon(aes(ymin = IIE_lwr, ymax = IIE_upr), alpha = 0.08, color = NA) +
	      geom_line(linewidth = 0.85) +
	      geom_point(size = 1.2) +
	      facet_wrap(~ contrast, scales = "free_y") +
	      labs(
	        x = "Residual correlation for each mediator varied one at a time",
	        y = "Mediator-specific Indirect Effect, posterior mean with 95% CrI",
	        color = "Mediator",
	        fill = "Mediator"
	      ) +
	      theme_bw(base_size = 12)
	  })

		  output$rho0_validation <- renderPrint({
	    s <- sensitivity()
    list(
      common_rho0_matches_main = s$rho0_common$passed,
      one_at_a_time_rho0_matches_main = s$rho0_one_at_a_time$passed
    )
  })

  sensitivity_selected <- reactive({
    req(identical(input$sensitivity_framework, "frequentist"))
    req(input$sensitivity_scenario, input$sensitivity_active_rho)

    if (identical(input$sensitivity_scenario, "common")) {
      sensitivity()$common$sensitivity_summary |>
        filter(valid_parameter, abs(rho - input$sensitivity_active_rho) < 1e-10) |>
        mutate(contrast = contrast_label_current(exposure)) |>
        filter_sensitivity_contrast() |>
        transmute(
          contrast,
          rho,
          adjusted_total_indirect = NIE_adjusted,
          adjusted_direct_effect = NDE_adjusted,
          adjusted_total_effect = NDE_adjusted + NIE_adjusted,
          PM_adjusted
        )
    } else {
      req(input$sensitivity_one_mediator)
      sensitivity()$one_at_a_time$sensitivity_summary |>
        filter(
          valid_parameter,
          varied_mediator == input$sensitivity_one_mediator,
          abs(rho - input$sensitivity_active_rho) < 1e-10
        ) |>
        mutate(contrast = contrast_label_current(exposure)) |>
        filter_sensitivity_contrast() |>
        transmute(
          contrast,
          varied_mediator,
          rho,
          adjusted_total_indirect = NIE_adjusted,
          adjusted_direct_effect = NDE_adjusted,
          adjusted_total_effect = NDE_adjusted + NIE_adjusted,
          PM_adjusted
        )
    }
  })

  output$sensitivity_selected <- renderDT({
    sensitivity_selected() |>
      present_frequentist_sensitivity_selected() |>
      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "t"))
  }, server = FALSE)

  output$frequentist_sensitivity_warnings <- renderUI({
    selected <- sensitivity_selected()
    warnings <- list()

    if (nrow(selected) == 0) {
      warnings <- c(warnings, list(tags$p("The selected rho is not feasible for the current frequentist residual covariance structure.")))
    }
    if (pm_warning_needed_frequentist(selected)) {
      warnings <- c(
        warnings,
        list(tags$p("Proportion mediated is difficult to interpret at this sensitivity setting because the total effect is near zero or the direct and indirect effects have opposing signs."))
      )
    }
    if (length(warnings) == 0) {
      return(NULL)
    }
    tags$div(style = "color:#8a5a00;", warnings)
  })

  output$sensitivity_curve <- renderPlot({
    req(identical(input$sensitivity_framework, "frequentist"))
    req(input$sensitivity_scenario)

    if (identical(input$sensitivity_scenario, "common")) {
      sensitivity()$common$sensitivity_summary |>
        filter(valid_parameter) |>
        mutate(contrast = contrast_label_current(exposure)) |>
        filter_sensitivity_contrast() |>
        ggplot(aes(x = rho, y = NIE_adjusted, color = contrast)) +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
        geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
        geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45") +
        geom_line(linewidth = 0.9) +
        geom_point(size = 1.4) +
        labs(x = "Residual correlation (rho)", y = "Total Indirect Effect", color = "Contrast") +
        theme_bw(base_size = 12)
    } else {
      req(input$sensitivity_one_mediator)
      sensitivity()$one_at_a_time$sensitivity_summary |>
        filter(valid_parameter, varied_mediator == input$sensitivity_one_mediator) |>
        mutate(contrast = contrast_label_current(exposure)) |>
        filter_sensitivity_contrast() |>
        ggplot(aes(x = rho, y = NIE_adjusted, color = contrast)) +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
        geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
        geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45") +
        geom_line(linewidth = 0.9) +
        geom_point(size = 1.4) +
        labs(x = paste("Residual correlation for", mediator_label(input$sensitivity_one_mediator)), y = "Total Indirect Effect", color = "Contrast") +
        theme_bw(base_size = 12)
    }
  })

  output$sensitivity_target_curve <- renderPlot({
    req(identical(input$sensitivity_framework, "frequentist"))
    req(identical(input$sensitivity_scenario, "one_at_a_time"))
    req(input$sensitivity_one_mediator)

    target <- input$sensitivity_one_mediator
    curve <- sensitivity()$one_at_a_time$sensitivity_path |>
      filter(
        valid_parameter,
        varied_mediator == target,
        mediator == target
      ) |>
      mutate(contrast = contrast_label_current(exposure)) |>
      filter_sensitivity_contrast()
    validate(need(nrow(curve) > 0, "No target mediator sensitivity curve rows match the current filters."))

    ggplot(curve, aes(x = rho, y = indirect_adjusted, color = contrast)) +
      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
      geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45") +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.4) +
      labs(
        x = paste("Residual correlation for", mediator_label(target)),
        y = paste(mediator_label(target), "Indirect Effect"),
        color = "Contrast"
      ) +
      theme_bw(base_size = 12)
  })

	  output$sensitivity_path_curve <- renderPlot({
	    req(identical(input$sensitivity_framework, "frequentist"))
	    req(input$sensitivity_scenario)

    if (identical(input$sensitivity_scenario, "common")) {
      sensitivity()$common$sensitivity_path |>
        filter(valid_parameter) |>
        mutate(contrast = contrast_label_current(exposure), mediator = mediator_label(mediator)) |>
        filter_sensitivity_contrast() |>
        ggplot(aes(x = rho, y = indirect_adjusted, color = mediator)) +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
        geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
        geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45") +
        geom_line(linewidth = 0.8) +
        facet_wrap(~ contrast, scales = "free_y") +
        labs(x = "Residual correlation (rho)", y = "Indirect Effect", color = "Mediator") +
        theme_bw(base_size = 12)
    } else {
      req(input$sensitivity_one_mediator)
      sensitivity()$one_at_a_time$sensitivity_path |>
        filter(valid_parameter, varied_mediator == input$sensitivity_one_mediator) |>
        mutate(
          contrast = contrast_label_current(exposure),
          mediator_id = mediator,
          mediator = mediator_label(mediator),
          role = if_else(mediator_id == input$sensitivity_one_mediator, "Varied mediator", "Other mediators")
        ) |>
        filter_sensitivity_contrast() |>
        ggplot(aes(x = rho, y = indirect_adjusted, color = mediator)) +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
        geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
        geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45") +
        geom_line(aes(linetype = role, alpha = role), linewidth = 0.8) +
        scale_alpha_manual(values = c("Varied mediator" = 1, "Other mediators" = 0.55), guide = "none") +
        facet_wrap(~ contrast, scales = "free_y") +
	        labs(x = paste("Residual correlation for", mediator_label(input$sensitivity_one_mediator)), y = "Indirect Effect", color = "Mediator", linetype = "Role") +
	        theme_bw(base_size = 12)
	    }
	  })

	  output$sensitivity_all_mediator_comparison <- renderPlot({
	    req(identical(input$sensitivity_framework, "frequentist"))
	    req(identical(input$sensitivity_scenario, "one_at_a_time"))

	    curve <- sensitivity()$one_at_a_time$sensitivity_path |>
	      filter(valid_parameter, mediator == varied_mediator) |>
	      mutate(
	        contrast = contrast_label_current(exposure),
	        mediator = mediator_label(mediator)
	      ) |>
	      filter_sensitivity_contrast()
	    validate(need(nrow(curve) > 0, "No mediator comparison rows match the current contrast filter."))

	    ggplot(curve, aes(x = rho, y = indirect_adjusted, color = mediator)) +
	      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
	      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
	      geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45") +
	      geom_line(linewidth = 0.85) +
	      geom_point(size = 1.2) +
	      facet_wrap(~ contrast, scales = "free_y") +
	      labs(
	        x = "Residual correlation for each mediator varied one at a time",
	        y = "Mediator-specific Indirect Effect",
	        color = "Mediator"
	      ) +
	      theme_bw(base_size = 12)
	  })

	  output$sensitivity_tipping <- renderDT({
    req(identical(input$sensitivity_framework, "frequentist"))
    req(input$sensitivity_scenario)

    tipping_tbl <- if (identical(input$sensitivity_scenario, "common")) {
      sensitivity()$tipping_common$total |>
        mutate(contrast = contrast_label_current(exposure)) |>
        filter_sensitivity_contrast() |>
        select(contrast, rho_zero_crossing, crossing_in_grid)
    } else {
      req(input$sensitivity_one_mediator)
      sensitivity()$tipping_one_at_a_time$total |>
        filter(varied_mediator == input$sensitivity_one_mediator) |>
        mutate(contrast = contrast_label_current(exposure), varied_mediator = mediator_label(varied_mediator)) |>
        filter_sensitivity_contrast() |>
        select(contrast, varied_mediator, rho_zero_crossing, crossing_in_grid)
    }

    tipping_tbl |>
      rename(
        `rho zero crossing` = rho_zero_crossing,
        `crossing in grid` = crossing_in_grid
      ) |>
      round_numeric(4) |>
      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "t"))
  }, server = FALSE)
}

shinyApp(ui = ui, server = server)
