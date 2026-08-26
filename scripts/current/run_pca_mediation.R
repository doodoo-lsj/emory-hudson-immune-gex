source("R/pca_mediation_pipeline.R")

results <- run_pca_mediation_pipeline(
  rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
  coord_x = "imagecol",
  coord_y = "imagerow",
  bootstrap_B = 1000,
  bootstrap_seed = 123,
  make_plots = TRUE
)

invisible(results)
