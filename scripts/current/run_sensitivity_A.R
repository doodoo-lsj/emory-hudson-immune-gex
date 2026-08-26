suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

source("R/pca_mediation_pipeline.R")
source("R/sensitivity_v2.R")
source("R/sensitivity_A_residual_rho.R")

run_sensitivity_A <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    coord_x = "imagecol",
    coord_y = "imagerow",
    rho_grid = seq(-0.5, 0.5, by = 0.05),
    tolerance = 1e-8) {
  analysis_input <- prepare_current_sensitivity_data(
    rds_path = rds_path,
    coord_x = coord_x,
    coord_y = coord_y
  )

  sensitivity_A_results <- run_sensitivity_A_residual_rho(
    data = analysis_input$dat_full_allpc,
    rho_grid = rho_grid,
    tolerance = tolerance
  )

  cat("\nSensitivity A mathematical status\n")
  print(sensitivity_A_results$method_specification)

  cat("\nBaseline decomposition\n")
  print(sensitivity_A_results$baseline$summary |> select(contrast, TE, direct_effect = NDE, total_IIE = NIE, PM))

  cat("\nBaseline mediator residual covariance Sigma_M\n")
  print(sensitivity_A_results$baseline$Sigma_M_unbiased)

  cat("\nBaseline mediator residual correlation\n")
  print(sensitivity_A_results$baseline$residual_correlation)

  cat("\nrho = 0 validation\n")
  print(sensitivity_A_results$validation[c("passed")])

  cat("\nCommon-rho tipping points: total IIE\n")
  print(sensitivity_A_results$sensitivity_A_common_rho$tipping$total)

  cat("\nOne-at-a-time tipping points: total IIE\n")
  print(sensitivity_A_results$sensitivity_A_one_at_a_time$tipping$total)

  cat("\nAvailable plot objects:\n")
  cat("  sensitivity_A_results$plots$total_common\n")
  cat("  sensitivity_A_results$plots$mediator_specific_common\n")
  cat("  sensitivity_A_results$plots$total_one_at_a_time\n")
  cat("  sensitivity_A_results$plots$mediator_specific_one_at_a_time\n")

  cat("\nPrimary result object: sensitivity_A_results\n")
  sensitivity_A_results
}

if (sys.nframe() == 0) {
  sensitivity_A_results <- run_sensitivity_A()
}
