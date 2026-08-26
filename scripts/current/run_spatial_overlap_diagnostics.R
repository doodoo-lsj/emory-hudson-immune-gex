suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source("R/pca_mediation_pipeline.R")
source("R/spatial_overlap_diagnostics.R")

run_spatial_overlap <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    coord_x = "imagecol",
    coord_y = "imagerow",
    spline_k = 50,
    spline_basis = "tp",
    grid_n = 120) {
  spatial_overlap_diagnostics <- run_spatial_overlap_diagnostics(
    rds_path = rds_path,
    coord_x = coord_x,
    coord_y = coord_y,
    spline_k = spline_k,
    spline_basis = spline_basis,
    grid_n = grid_n
  )

  cat("\nSpline specification\n")
  print(spatial_overlap_diagnostics$config[c("spline_specification", "method")])

  cat("\nAnalysis data predictor scales\n")
  print(spatial_overlap_diagnostics$data$predictor_scale_summary)

  cat("\nCoordinate standardization used by the fixed spline specification\n")
  print(spatial_overlap_diagnostics$data$coordinate_scale_info)

  cat("\nGAM fit summary\n")
  print(spatial_overlap_diagnostics$gam_fit_summary)

  cat("\nSmooth EDF and tests\n")
  print(spatial_overlap_diagnostics$smooth_edf)

  cat("\nk-index / basis adequacy diagnostics\n")
  print(spatial_overlap_diagnostics$k_diagnostics)

  cat("\nConcurvity summary, full = TRUE\n")
  print(spatial_overlap_diagnostics$concurvity$tidy_full)

  cat("\nConcurvity summary, full = FALSE\n")
  print(spatial_overlap_diagnostics$concurvity$tidy_pairwise)

  cat("\nExposure spatial predictability: X dummy ~ s(x_std, y_std)\n")
  print(spatial_overlap_diagnostics$exposure_spatial_predictability$summary)

  cat("\nPredicted probability distribution by observed exposure category\n")
  print(spatial_overlap_diagnostics$exposure_spatial_predictability$probability_distribution)

  cat("\nMediator spatial smooth contribution by exposure category\n")
  print(spatial_overlap_diagnostics$smooth_by_exposure$group_summary)

  cat("\nMediator spatial smooth contribution group mean differences\n")
  print(spatial_overlap_diagnostics$smooth_by_exposure$group_mean_differences)

  cat("\nLinear vs spline alpha stability\n")
  print(spatial_overlap_diagnostics$alpha_stability)

  cat("\nPlot objects are available at:\n")
  cat("  spatial_overlap_diagnostics$plots$exposure_map\n")
  cat("  spatial_overlap_diagnostics$plots$adjacent_probability_surface\n")
  cat("  spatial_overlap_diagnostics$plots$far_probability_surface\n")

  cat("\nPrimary result object: spatial_overlap_diagnostics\n")
  spatial_overlap_diagnostics
}

## Rscript scripts/current/run_spatial_overlap_diagnostics.R executes the diagnostics.
## source("scripts/current/run_spatial_overlap_diagnostics.R") only defines run_spatial_overlap().
if (sys.nframe() == 0) {
  spatial_overlap_diagnostics <- run_spatial_overlap()
}
