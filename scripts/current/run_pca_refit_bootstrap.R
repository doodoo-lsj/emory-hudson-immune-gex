suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source("R/pca_mediation_pipeline.R")
source("R/pca_refit_bootstrap.R")

run_pca_refit_bootstrap <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    B = 1000,
    seed = 123,
    n_cores = 1,
    min_loading_correlation = 0.70,
    run_fixed_pca_bootstrap = TRUE) {
  pca_refit_bootstrap_results <- run_pca_refit_bootstrap_analysis(
    rds_path = rds_path,
    B = B,
    seed = seed,
    n_cores = n_cores,
    min_loading_correlation = min_loading_correlation,
    run_fixed_pca_bootstrap = run_fixed_pca_bootstrap
  )

  cat("\nPCA-refit bootstrap config\n")
  print(pca_refit_bootstrap_results$config)

  cat("\nAlignment helper sanity check\n")
  print(pca_refit_bootstrap_results$validation$alignment_sanity_check$passed)

  cat("\nFailure summary\n")
  print(pca_refit_bootstrap_results$failure_summary |> select(-failure_details))

  cat("\nAlignment summary\n")
  print(
    pca_refit_bootstrap_results$alignment_summary |>
      summarise(
        min_matched_abs_correlation_min = min(min_matched_abs_correlation, na.rm = TRUE),
        min_matched_abs_correlation_median = median(min_matched_abs_correlation, na.rm = TRUE),
        mean_matched_abs_correlation_mean = mean(mean_matched_abs_correlation, na.rm = TRUE),
        alignment_warning_B = sum(alignment_warning, na.rm = TRUE)
      )
  )

  cat("\nExplained variance stability for matched bootstrap components\n")
  print(pca_refit_bootstrap_results$explained_variance_stability)

  cat("\nLoading alignment stability\n")
  print(pca_refit_bootstrap_results$loading_alignment_stability)

  if (!is.null(pca_refit_bootstrap_results$fixed_vs_refit_comparison)) {
    cat("\nFixed-PCA vs PCA-refit bootstrap comparison, all successful refit replicates\n")
    print(pca_refit_bootstrap_results$fixed_vs_refit_comparison)

    cat("\nFixed-PCA vs PCA-refit bootstrap comparison, high-quality aligned refit replicates only\n")
    print(pca_refit_bootstrap_results$fixed_vs_refit_comparison_high_quality)
  }

  cat("\nTiming\n")
  print(pca_refit_bootstrap_results$timing)

  cat("\nPrimary result object: pca_refit_bootstrap_results\n")
  pca_refit_bootstrap_results
}

## Rscript scripts/current/run_pca_refit_bootstrap.R executes B = 1000 by default.
## source("scripts/current/run_pca_refit_bootstrap.R") only defines run_pca_refit_bootstrap().
if (sys.nframe() == 0) {
  pca_refit_bootstrap_results <- run_pca_refit_bootstrap()
}
