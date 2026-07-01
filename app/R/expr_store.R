# expr_store.R — multi-dataset HDF5 backing store for per-gene expression reads.
#
# Each dataset has its own app/data/<key>/expression.h5 holding "expr"
# (n_genes x n_cells, log-norm, one gene-row per chunk), "genes" (row order) and
# "cells" (column order == that dataset's umap_metadata row order). We read each
# dataset's gene index once at startup and pull a single gene-row on demand via
# an rhdf5 hyperslab, so no full matrix ever lives in memory.
#
# rhdf5 (not hdf5r) is used because it bundles its own HDF5 library — self-
# contained for both local runs and shinyapps.io deployment. It ships in the
# base R 4.5.3 library here, so no Seurat/renv is needed at runtime.

.expr_env <- new.env(parent = emptyenv())

# init_expr_stores(keys) — open every dataset's gene index once.
init_expr_stores <- function(keys) {
  for (k in keys) {
    path <- file.path(app_data_dir, k, "expression.h5")
    if (!file.exists(path)) next
    genes <- as.character(rhdf5::h5read(path, "genes"))
    .expr_env[[k]] <- list(
      path  = path,
      index = stats::setNames(seq_along(genes), genes)
    )
  }
  invisible(TRUE)
}

# read_gene(dataset, gene) -> numeric vector length n_cells aligned to that
# dataset's "cells" / umap_metadata row order, or NULL if absent / unreadable.
#
# Orientation note (as in kpc_tme): h5py writes "expr" as (genes x cells) but
# rhdf5 reads it back transposed to (cells x genes), so the gene is the SECOND
# index: list(NULL, idx) selects all cells for one gene, aligned to "cells".
read_gene <- function(dataset, gene) {
  st <- .expr_env[[dataset]]
  if (is.null(st)) return(NULL)
  idx <- st$index[gene]
  if (is.na(idx)) return(NULL)
  idx <- as.integer(idx)
  tryCatch(
    as.numeric(rhdf5::h5read(st$path, "expr", index = list(NULL, idx))),
    error = function(e) NULL
  )
}

# read_genes(dataset, genes) -> named list gene -> vector (skips absent genes).
# Used by the marker heatmap / dot plot which aggregate several genes at once.
read_genes <- function(dataset, genes) {
  out <- lapply(genes, function(g) read_gene(dataset, g))
  names(out) <- genes
  out[!vapply(out, is.null, logical(1))]
}
