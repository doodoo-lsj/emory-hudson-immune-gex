suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source("R/pca_mediation_pipeline.R")
source("R/spatial_sensitivity.R")

run_spatial_sensitivity <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    coord_x = "imagecol",
    coord_y = "imagerow",
    spline_k = 50,
    spline_basis = "tp",
    include = c("linear", "quadratic", "spline"),
    compute_moran = TRUE,
    moran_k = 6) {
  spatial_comparison_results <- run_spatial_sensitivity_analysis(
    rds_path = rds_path,
    coord_x = coord_x,
    coord_y = coord_y,
    spline_k = spline_k,
    spline_basis = spline_basis,
    include = include,
    compute_moran = compute_moran,
    moran_k = moran_k
  )

  cat("\nAnalysis data predictor scales\n")
  print(spatial_comparison_results$data$predictor_scale_summary)

  cat("\nCoordinate standardization used for quadratic and spline spatial terms\n")
  print(spatial_comparison_results$data$coordinate_scale_info)

  cat("\nSpatial specifications\n")
  print(spatial_comparison_results$specifications)

  cat("\nMediator residual pairwise correlations\n")
  print(
    spatial_comparison_results$residual_pairwise |>
      select(spatial_specification, mediator_1, mediator_2, residual_correlation)
  )

  cat("\nAlpha comparison: X -> mediator paths\n")
  print(
    spatial_comparison_results$alpha_comparison |>
      select(
        spatial_specification,
        contrast,
        mediator,
        alpha,
        SE,
        alpha_absolute_change,
        alpha_relative_change
      )
  )

  cat("\nBeta comparison: mediator -> outcome paths\n")
  print(
    spatial_comparison_results$beta_comparison |>
      select(
        spatial_specification,
        mediator,
        beta,
        SE,
        beta_absolute_change,
        beta_relative_change
      )
  )

  cat("\nMediator-specific indirect component comparison\n")
  print(
    spatial_comparison_results$indirect_comparison |>
      select(
        spatial_specification,
        contrast,
        mediator,
        alpha,
        beta,
        indirect_component,
        indirect_component_absolute_change,
        indirect_component_relative_change
      )
  )

  cat("\nTotal mediation decomposition comparison\n")
  print(
    spatial_comparison_results$decomposition_comparison |>
      select(
        spatial_specification,
        contrast,
        TE,
        direct_effect,
        total_NIE,
        PM,
        total_NIE_absolute_change,
        total_NIE_relative_change
      )
  )

  if (!is.null(spatial_comparison_results$spatial_diagnostics)) {
    cat("\nOptional k-nearest-neighbor Moran's I residual diagnostic\n")
    print(spatial_comparison_results$spatial_diagnostics)
  }

  cat("\nPrimary result object: spatial_comparison_results\n")
  spatial_comparison_results
}

## Rscript scripts/current/run_spatial_sensitivity.R executes the analysis.
## source("scripts/current/run_spatial_sensitivity.R") only defines run_spatial_sensitivity().
if (sys.nframe() == 0) {
  spatial_comparison_results <- run_spatial_sensitivity()
}
