suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

source("R/pca_mediation_pipeline.R")
source("R/sensitivity_A_residual_rho.R")
source("R/sensitivity_B_latent_U.R")

run_sensitivity_B <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    coord_x = "imagecol",
    coord_y = "imagerow",
    p_U_grid = 0.3,
    OR_UX_grid = c(1, 1.5, 2, 3),
    s_M_grid = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5),
    s_Y_grid = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5),
    scenarios = c("B1_XY", "B2_XM", "B3_MY", "B4_global"),
    maxit = 1000,
    reltol = 1e-8,
    null_tolerance = 1e-5) {
  analysis_input <- prepare_current_sensitivity_data(
    rds_path = rds_path,
    coord_x = coord_x,
    coord_y = coord_y
  )

  sensitivity_B_results <- run_sensitivity_B_latent_U(
    data = analysis_input$dat_full_allpc,
    p_U_grid = p_U_grid,
    OR_UX_grid = OR_UX_grid,
    s_M_grid = s_M_grid,
    s_Y_grid = s_Y_grid,
    scenarios = scenarios,
    maxit = maxit,
    reltol = reltol,
    null_tolerance = null_tolerance
  )

  cat("\nSensitivity B observed-data likelihood\n")
  print(sensitivity_B_results$likelihood_specification)

  cat("\nBaseline decomposition\n")
  print(sensitivity_B_results$baseline$summary |> select(contrast, TE, direct_effect = NDE, total_IIE = NIE, PM))

  cat("\nNull validation\n")
  print(sensitivity_B_results$null_validation[c("passed", "tolerance")])

  cat("\nScenario-specific null mechanism validation\n")
  print(sensitivity_B_results$scenario_null_validation)

  cat("\nB3/B4 nesting validation\n")
  print(sensitivity_B_results$nesting_validation)

  cat("\nTE decomposition identity validation\n")
  print(sensitivity_B_results$decomposition_validation)

  cat("\nPath behavior validation\n")
  print(sensitivity_B_results$path_behavior_validation)

  cat("\nRepresentative B3 path changes\n")
  print(
    sensitivity_B_results$representative_B3 |>
      select(
        s_M1, s_Y, contrast, mediator, alpha, beta, IIE,
        baseline_beta, beta_change, baseline_IIE, IIE_change,
        max_abs_beta_change_scenario, max_abs_IIE_change_scenario
      )
  )

  cat("\nRuntime by scenario\n")
  print(sensitivity_B_results$timing)

  cat("\nConvergence summaries\n")
  for (scenario in c("sensitivity_B1_XY", "sensitivity_B2_XM", "sensitivity_B3_MY", "sensitivity_B4_global")) {
    if (!is.null(sensitivity_B_results[[scenario]])) {
      cat("\n", scenario, "\n", sep = "")
      convergence_tbl <- sensitivity_B_results[[scenario]]$convergence
      if (!"failure_reason" %in% names(convergence_tbl)) {
        convergence_tbl$failure_reason <- NA_character_
      }
      print(
        convergence_tbl |>
          summarise(
            grid_points = dplyr::n(),
            joint_nonzero_convergence = sum(joint_convergence != 0, na.rm = TRUE),
            total_nonzero_convergence = sum(total_convergence != 0, na.rm = TRUE),
            failures = sum(!is.na(.data$failure_reason), na.rm = TRUE)
          )
      )
    }
  }

  cat("\nPrimary result object: sensitivity_B_results\n")
  sensitivity_B_results
}

## Example full-ish run:
## sensitivity_B_results <- run_sensitivity_B(
##   p_U_grid = 0.3,
##   OR_UX_grid = c(1, 1.5, 2, 3),
##   s_M_grid = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5),
##   s_Y_grid = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5)
## )

if (sys.nframe() == 0) {
  sensitivity_B_results <- run_sensitivity_B()
}
