suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

source("R/sensitivity_v2.R")

rds_path <- "data/raw/pt16_emory_GEX_immune_FULL_v2.rds"
coord_x <- "imagecol"
coord_y <- "imagerow"

exposures <- c("X_adjacent", "X_far")
mediators <- c("PC1_R", "PC2_R", "PC3")
covariates <- c("x_coord", "y_coord")

# Sensitivity parameter:
# rho_j = Corr(error_Mj, error_Y), where error_Mj is from
# M_j ~ X_adjacent + X_far + x_coord + y_coord and error_Y is from
# Y ~ X_adjacent + X_far + PC1_R + PC2_R + PC3 + x_coord + y_coord.
rho_grid <- seq(-0.5, 0.5, by = 0.05)

analysis_input <- prepare_v2_mediation_data(
  rds_path = rds_path,
  coord_x = coord_x,
  coord_y = coord_y
)

dat_full_allpc <- analysis_input$dat_full_allpc

main_decomposition <- decompose_linear_multix(
  data = dat_full_allpc,
  exposures = exposures,
  outcome = "Y",
  mediators = mediators,
  covariates = covariates
)

sensitivity_common <- residual_corr_sensitivity(
  data = dat_full_allpc,
  exposures = exposures,
  outcome = "Y",
  mediators = mediators,
  covariates = covariates,
  rho_grid = rho_grid,
  mode = "common"
)

sensitivity_one_at_a_time <- residual_corr_sensitivity(
  data = dat_full_allpc,
  exposures = exposures,
  outcome = "Y",
  mediators = mediators,
  covariates = covariates,
  rho_grid = rho_grid,
  mode = "one_at_a_time"
)

rho0_check_common <- validate_rho0_matches_observed(sensitivity_common)
rho0_check_one_at_a_time <- validate_rho0_matches_observed(sensitivity_one_at_a_time)

if (!isTRUE(rho0_check_common$passed) || !isTRUE(rho0_check_one_at_a_time$passed)) {
  stop("rho = 0 sensitivity estimates do not match the main mediation estimates.")
}

add_contrast_label <- function(df) {
  df |>
    mutate(
      contrast = case_when(
        exposure == "X_adjacent" ~ "adjacent vs inside",
        exposure == "X_far" ~ "far vs inside",
        TRUE ~ exposure
      )
    )
}

original_summary <- sensitivity_common$observed_summary |>
  add_contrast_label() |>
  select(contrast, exposure, TE, direct_effect = NDE, total_indirect_effect = NIE, PM)

original_path <- sensitivity_common$observed_path |>
  add_contrast_label() |>
  select(
    contrast,
    exposure,
    mediator,
    alpha = alpha_X_to_M,
    beta = beta_M_to_Y,
    indirect_component
  )

sensitivity_total_common <- sensitivity_common$sensitivity_summary |>
  add_contrast_label() |>
  select(
    scenario,
    contrast,
    exposure,
    rho,
    starts_with("rho_"),
    adjusted_total_indirect = NIE_adjusted,
    adjusted_direct_effect = NDE_adjusted,
    PM_adjusted,
    feasibility_ratio,
    valid_parameter
  )

sensitivity_path_common <- sensitivity_common$sensitivity_path |>
  add_contrast_label() |>
  select(
    scenario,
    contrast,
    exposure,
    mediator,
    rho,
    starts_with("rho_"),
    alpha = alpha_X_to_M,
    beta_observed,
    adjusted_beta = beta_adjusted,
    adjusted_indirect_component = indirect_adjusted,
    feasibility_ratio,
    valid_parameter
  )

sensitivity_total_one_at_a_time <- sensitivity_one_at_a_time$sensitivity_summary |>
  add_contrast_label() |>
  select(
    scenario,
    varied_mediator,
    contrast,
    exposure,
    rho,
    starts_with("rho_"),
    adjusted_total_indirect = NIE_adjusted,
    adjusted_direct_effect = NDE_adjusted,
    PM_adjusted,
    feasibility_ratio,
    valid_parameter
  )

sensitivity_path_one_at_a_time <- sensitivity_one_at_a_time$sensitivity_path |>
  add_contrast_label() |>
  select(
    scenario,
    varied_mediator,
    contrast,
    exposure,
    mediator,
    rho,
    starts_with("rho_"),
    alpha = alpha_X_to_M,
    beta_observed,
    adjusted_beta = beta_adjusted,
    adjusted_indirect_component = indirect_adjusted,
    feasibility_ratio,
    valid_parameter
  )

tipping_common <- summarize_rho_tipping_points(sensitivity_common)
tipping_one_at_a_time <- summarize_rho_tipping_points(sensitivity_one_at_a_time)

sensitivity_plots <- build_rho_sensitivity_plots(
  common_result = sensitivity_common,
  one_at_a_time_result = sensitivity_one_at_a_time
)

validity_summary <- list(
  rho0_common_passed = rho0_check_common$passed,
  rho0_one_at_a_time_passed = rho0_check_one_at_a_time$passed,
  common_valid_rho_range = range(
    sensitivity_total_common$rho[sensitivity_total_common$valid_parameter],
    na.rm = TRUE
  ),
  common_invalid_rows = sum(!sensitivity_total_common$valid_parameter),
  one_at_a_time_valid_rho_range = range(
    sensitivity_total_one_at_a_time$rho[sensitivity_total_one_at_a_time$valid_parameter],
    na.rm = TRUE
  ),
  one_at_a_time_invalid_rows = sum(!sensitivity_total_one_at_a_time$valid_parameter)
)

cat("\nOriginal coordinate-adjusted mediation estimates\n")
print(original_summary)

cat("\nOriginal mediator-specific path estimates\n")
print(original_path)

cat("\nrho = 0 validation\n")
print(validity_summary[c("rho0_common_passed", "rho0_one_at_a_time_passed")])

cat("\nFeasible rho ranges for the current grid\n")
print(validity_summary[c(
  "common_valid_rho_range",
  "common_invalid_rows",
  "one_at_a_time_valid_rho_range",
  "one_at_a_time_invalid_rows"
)])

cat("\nCommon-rho total NIE tipping points\n")
print(tipping_common$total)

cat("\nOne-at-a-time total NIE tipping points\n")
print(tipping_one_at_a_time$total)

cat("\nAvailable plot objects: sensitivity_plots$total_common, ",
    "sensitivity_plots$mediator_specific_common, ",
    "sensitivity_plots$total_one_at_a_time, ",
    "sensitivity_plots$mediator_specific_one_at_a_time\n",
    sep = "")

invisible(list(
  dat_full_allpc = dat_full_allpc,
  main_decomposition = main_decomposition,
  original_summary = original_summary,
  original_path = original_path,
  sensitivity_common = sensitivity_common,
  sensitivity_one_at_a_time = sensitivity_one_at_a_time,
  sensitivity_total_common = sensitivity_total_common,
  sensitivity_path_common = sensitivity_path_common,
  sensitivity_total_one_at_a_time = sensitivity_total_one_at_a_time,
  sensitivity_path_one_at_a_time = sensitivity_path_one_at_a_time,
  tipping_common = tipping_common,
  tipping_one_at_a_time = tipping_one_at_a_time,
  sensitivity_plots = sensitivity_plots,
  validity_summary = validity_summary
))
