suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source("R/pca_mediation_pipeline.R")
source("R/joint_mediator_frequentist.R")

run_joint_mediator_frequentist <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    coord_x = "imagecol",
    coord_y = "imagerow",
    standardize_covariates = FALSE,
    continuous_covariates = c("x_coord", "y_coord"),
    run_standardization_check = TRUE,
    tolerance = 1e-8) {
  joint_mediator_results <- run_joint_mediator_frequentist_analysis(
    rds_path = rds_path,
    coord_x = coord_x,
    coord_y = coord_y,
    standardize_covariates = standardize_covariates,
    continuous_covariates = continuous_covariates,
    run_standardization_check = run_standardization_check,
    tolerance = tolerance
  )

  cat("\nAnalysis data predictor scales\n")
  print(joint_mediator_results$data$predictor_scale_summary)

  if (isTRUE(standardize_covariates)) {
    cat("\nContinuous covariates standardized for fitting\n")
    print(joint_mediator_results$data$scale_info)
    cat("\nCoefficient tables and mediation quantities are returned on the original covariate scale.\n")
  }

  cat("\nSeparate OLS vs joint multivariate mediator coefficient comparison\n")
  print(
    joint_mediator_results$coefficient_comparison |>
      filter(term %in% c("X_adjacent", "X_far", "x_coord", "y_coord"))
  )

  cat("\nCoefficient equivalence validation\n")
  print(joint_mediator_results$validation)

  cat("\nJoint mediator residual covariance matrix: Sigma_M_hat = E'E / n\n")
  print(joint_mediator_results$residual_covariance)

  cat("\nJoint mediator residual correlation matrix\n")
  print(joint_mediator_results$residual_correlation)

  cat("\nResidual pairwise summary\n")
  print(joint_mediator_results$residual_pairwise_summary)

  cat("\nAlpha cross-equation coefficient covariance: X_adjacent\n")
  print(joint_mediator_results$alpha_cross_equation_covariance$X_adjacent$covariance)

  cat("\nAlpha cross-equation coefficient covariance: X_far\n")
  print(joint_mediator_results$alpha_cross_equation_covariance$X_far$covariance)

  cat("\nMediator-specific indirect comparison\n")
  print(
    joint_mediator_results$mediation_comparison$path |>
      select(
        contrast,
        mediator,
        original_indirect,
        joint_model_indirect,
        indirect_difference
      )
  )

  cat("\nTotal NIE comparison\n")
  print(
    joint_mediator_results$mediation_comparison$summary |>
      select(
        contrast,
        original_NIE,
        joint_model_NIE,
        NIE_difference
      )
  )

  if (!is.null(joint_mediator_results$standardization_check)) {
    cat("\nStandardized vs unstandardized mediation point-estimate check\n")
    print(joint_mediator_results$standardization_check$summary)
    cat("\nPassed: ", joint_mediator_results$standardization_check$passed, "\n", sep = "")
  }

  cat("\nPrimary result object: joint_mediator_results\n")
  joint_mediator_results
}

## Rscript scripts/current/run_joint_mediator_frequentist.R executes the analysis.
## source("scripts/current/run_joint_mediator_frequentist.R") only defines run_joint_mediator_frequentist().
if (sys.nframe() == 0) {
  joint_mediator_results <- run_joint_mediator_frequentist()
}
