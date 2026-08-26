suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source("R/pca_mediation_pipeline.R")
source("R/pca_refit_bootstrap.R")
source("R/pca_axis_stability_diagnostics.R")

run_pca_axis_stability <- function(
    results_path = "pca_refit_bootstrap.RData",
    object_name = "pca_refit_bootstrap_results",
    recompute_full_similarity_matrices = FALSE,
    bootstrap_ids_for_recompute = NULL,
    n_cores = 1) {
  pca_axis_stability_diagnostics <- run_pca_axis_stability_diagnostics(
    results_path = results_path,
    object_name = object_name,
    recompute_full_similarity_matrices = recompute_full_similarity_matrices,
    bootstrap_ids_for_recompute = bootstrap_ids_for_recompute,
    n_cores = n_cores
  )

  cat("\nPCA axis stability config\n")
  print(pca_axis_stability_diagnostics$config)

  cat("\nAlignment method used by current PCA-refit bootstrap\n")
  print(pca_axis_stability_diagnostics$alignment_method)

  cat("\nAligned absolute loading correlation distribution by reference component\n")
  print(pca_axis_stability_diagnostics$axis_similarity_distribution)

  cat("\nPC permutation / swap summary\n")
  print(pca_axis_stability_diagnostics$permutation_summary$total)

  cat("\nPC permutation pattern counts\n")
  print(pca_axis_stability_diagnostics$permutation_summary$pattern_counts)

  cat("\nPer-axis matched bootstrap component counts\n")
  print(pca_axis_stability_diagnostics$permutation_summary$per_axis)

  cat("\nSign flip summary\n")
  print(pca_axis_stability_diagnostics$sign_flip_summary$total)

  cat("\nSign flips by reference component\n")
  print(pca_axis_stability_diagnostics$sign_flip_summary$per_axis)

  cat("\nExplained-variance component stability\n")
  print(pca_axis_stability_diagnostics$explained_variance_ordering$component_summary)

  cat("\nExplained-variance gap stability\n")
  print(pca_axis_stability_diagnostics$explained_variance_ordering$gap_summary)

  cat("\nAligned explained-variance ordering summary\n")
  print(pca_axis_stability_diagnostics$explained_variance_ordering$ordering_summary)

  cat("\nPC1-PC3 3D subspace stability by principal-angle dimension\n")
  print(pca_axis_stability_diagnostics$subspace_stability$per_dimension)

  cat("\nPC1-PC3 3D subspace replicate-level summary\n")
  print(pca_axis_stability_diagnostics$subspace_stability$replicate_summary)

  cat("\nIndividual-axis vs 3D-subspace summary\n")
  print(pca_axis_stability_diagnostics$axis_vs_subspace_summary)

  cat("\nFull 3x3 similarity matrix status\n")
  print(pca_axis_stability_diagnostics$full_similarity_note)

  if (!is.null(pca_axis_stability_diagnostics$full_similarity_matrices)) {
    cat("\nFull 3x3 loading correlation/cosine matrices\n")
    print(pca_axis_stability_diagnostics$full_similarity_matrices)
  }

  cat("\nPrimary result object: pca_axis_stability_diagnostics\n")
  pca_axis_stability_diagnostics
}

## Rscript scripts/current/run_pca_axis_stability_diagnostics.R runs diagnostics
## from the saved PCA-refit bootstrap object without recomputing full 3x3 matrices.
## source("scripts/current/run_pca_axis_stability_diagnostics.R") only defines run_pca_axis_stability().
if (sys.nframe() == 0) {
  pca_axis_stability_diagnostics <- run_pca_axis_stability()
}
