# Spatial Network Mediation Project

This repository contains the spatial transcriptomics mediation analysis for immune infiltration and tumor proximity.

## Current Entry Points

- Canonical frequentist pipeline: `scripts/current/pcaaxis_0722.R`
- Optional z-score mediator version: `scripts/optional/pcaaxis_0722_zscore.R`
- Current v2 sensitivity helper functions: `R/sensitivity_v2.R`
- Analysis decisions and audit notes: `docs/audit/`

## Current Data

- Current input data: `data/raw/pt16_emory_GEX_immune_FULL_v2.rds`
- Historical input data: `data/raw/pt16_emory_GEX_immune_FULL.rds`
- Metadata workbook: `data/raw/pt16_emory_GEX_immune_metadata.xlsx`

## Folder Map

- `R/`: reusable helper functions.
- `scripts/current/`: current runnable analysis scripts.
- `scripts/optional/`: optional analysis variants.
- `scripts/legacy/`: historical exploratory scripts; do not treat as current results.
- `data/raw/`: raw input data.
- `results/figures/current/`: current figures.
- `results/figures/legacy/`: legacy or provenance-uncertain figures.
- `results/enrichment/current/`: current ToppGene enrichment downloads.
- `results/tables/legacy/`: legacy result tables.
- `docs/meetings/`: meeting notes and meeting-preparation files.
- `docs/audit/`: audit reports, cleanup logs, and project decisions.
- `archive/before_data_update/`: pre-data-update enrichment archive.
- `archive/workspaces/`: archived R workspaces.

## Notes

The project folder is not currently a git repository. No files were deleted during cleanup.
