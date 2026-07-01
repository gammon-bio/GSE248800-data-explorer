#!/usr/bin/env Rscript
# 03_build_app_data.R — R half of the data prep.
#
# Consumes the per-dataset scaffolding written by 02_export_from_h5ad.py and
# builds the lightweight RDS artifacts the Shiny app loads per dataset:
#
#   app/data/<dataset>/umap_metadata.rds      (row order == expression.h5 cells)
#   app/data/<dataset>/gene_list.rds          (sorted, for the dropdown)
#   app/data/<dataset>/cell_type_palette.rds  (SHARED palette across datasets)
#   app/data/<dataset>/summary_stats.rds      (cell_type x condition counts)
#   app/data/<dataset>/de_results.rds         (named list: cell_type -> DE df)
#   app/data/<dataset>/marker_genes.rds       (one-vs-rest markers)
#
# expression.h5 and the DE/marker CSVs are written by the Python step; this
# script derives the .rds and verifies the heavy artifacts + invariants.
#
# Run:  Rscript data-raw/03_build_app_data.R   (uses base R 4.5.3, no Seurat)

suppressPackageStartupMessages({
  library(dplyr)
})

.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
REPO <- if (length(.file_arg) == 1L) {
  dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
} else normalizePath(getwd())

APP_DATA <- file.path(REPO, "app", "data")
DATASETS <- c("mouse_pancreatic", "mouse_colorectal", "human_pdac")

# Canonical lineage -> colour so cell types are coloured consistently across the
# three datasets (any lineage not listed gets a generated fallback colour).
LINEAGE_COLORS <- c(
  "Macrophage/Monocyte" = "#E41A1C",
  "Dendritic cell"      = "#FF7F00",
  "T/NK cell"           = "#377EB8",
  "B cell"              = "#4DAF4A",
  "Neutrophil"          = "#984EA3",
  "FAP/Fibroblast"      = "#A65628",
  "Endothelial"         = "#66C2A5",
  "Pericyte/Mural"      = "#999999",
  "Satellite cell"      = "#F781BF",
  "Adipocyte"           = "#FFC107",
  "Schwann/Glia"        = "#1B9E77",
  "Tenocyte"            = "#C49A6C"
)

read_scaffold_meta <- function(dir) {
  utils::read.csv(file.path(dir, "_scaffold", "umap_metadata.csv"),
                  stringsAsFactors = FALSE, colClasses = c(cell = "character"))
}

# ---- Pass 1: collect the union of cell types for a shared palette ----------
all_types <- character(0)
for (ds in DATASETS) {
  d <- file.path(APP_DATA, ds)
  if (!dir.exists(file.path(d, "_scaffold"))) next
  all_types <- union(all_types, unique(read_scaffold_meta(d)$cell_type))
}
all_types <- sort(all_types)
palette <- LINEAGE_COLORS[all_types]
missing <- all_types[is.na(palette)]
if (length(missing)) {
  fallback <- grDevices::hcl.colors(length(missing), palette = "Dark3")
  palette[is.na(palette)] <- fallback
}
names(palette) <- all_types
cat(sprintf("[03] shared palette: %d cell types\n", length(palette)))

# ---- Pass 2: per-dataset artifacts -----------------------------------------
for (ds in DATASETS) {
  d <- file.path(APP_DATA, ds)
  scaffold <- file.path(d, "_scaffold")
  if (!dir.exists(scaffold)) { cat(sprintf("[03] SKIP %s (no scaffold)\n", ds)); next }
  cat(sprintf("\n[03] ==== %s ====\n", ds))

  meta <- read_scaffold_meta(d)                       # DO NOT reorder rows
  saveRDS(meta, file.path(d, "umap_metadata.rds"))
  cat(sprintf("[03]   umap_metadata.rds: %d cells\n", nrow(meta)))

  genes <- readLines(file.path(scaffold, "genes.txt"))
  genes <- genes[nzchar(genes)]
  saveRDS(sort(unique(genes)), file.path(d, "gene_list.rds"))
  cat(sprintf("[03]   gene_list.rds: %d genes\n", length(genes)))

  saveRDS(palette, file.path(d, "cell_type_palette.rds"))

  summary_stats <- meta %>%
    dplyr::count(cell_type, condition, name = "n_cells") %>%
    dplyr::arrange(cell_type, condition) %>%
    as.data.frame()
  saveRDS(summary_stats, file.path(d, "summary_stats.rds"))
  cat(sprintf("[03]   summary_stats.rds: %d rows\n", nrow(summary_stats)))

  # de_results: named list keyed by cell_type (from de_*.csv 'cell_type' column)
  de_files <- list.files(d, pattern = "^de_.*\\.csv$", full.names = TRUE)
  de_results <- list()
  for (f in de_files) {
    df <- utils::read.csv(f, stringsAsFactors = FALSE)
    ct <- df$cell_type[1]
    df$cell_type <- NULL
    de_results[[ct]] <- df
  }
  saveRDS(de_results, file.path(d, "de_results.rds"))
  cat(sprintf("[03]   de_results.rds: %d cell types\n", length(de_results)))

  marker_path <- file.path(d, "marker_genes.csv")
  if (file.exists(marker_path)) {
    saveRDS(utils::read.csv(marker_path, stringsAsFactors = FALSE),
            file.path(d, "marker_genes.rds"))
    cat("[03]   marker_genes.rds written\n")
  }

  # ---- Verify heavy artifact + the cells<->metadata invariant --------------
  # Check against the "cells" dataset length directly (robust to rhdf5's
  # row/column-major dim reversal, which makes h5ls report cells x genes).
  h5 <- file.path(d, "expression.h5")
  if (file.exists(h5)) {
    n_cells_h5 <- tryCatch(length(rhdf5::h5read(h5, "cells")), error = function(e) NA_integer_)
    cat(sprintf("[03]   expression.h5 present (%d cells)\n", n_cells_h5))
    if (!is.na(n_cells_h5) && n_cells_h5 != nrow(meta)) {
      warning(sprintf("[03]   MISMATCH expr cells %d != meta rows %d", n_cells_h5, nrow(meta)))
    }
  } else {
    warning("[03]   MISSING expression.h5 for ", ds)
  }
}

cat("\n[03] === App data build complete ===\n")
