# global.R — app startup.
# Loads each dataset's lightweight RDS artifacts into the DS registry, opens the
# per-dataset HDF5 expression handles, and sources helpers + feature modules.
# The three datasets the app toggles between are small enough to hold their
# metadata/DE/markers in memory at once; only the expression matrices stay on
# disk (read one gene-row at a time via R/expr_store.R). Shiny auto-sources this
# file once at launch.

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(shinyWidgets)
  library(DT)
  library(plotly)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(scales)
  library(waiter)
  library(htmltools)
  library(rhdf5)
})

app_data_dir <- "data"

# Idle timeout (minutes). After this long with no user interaction the client
# closes its Shiny WebSocket so a forgotten open tab stops holding the server
# awake — letting a serverless host (Railway) sleep the container and not burn
# credits. Override per-deploy with the IDLE_MINUTES env var.
IDLE_MINUTES <- as.integer(Sys.getenv("IDLE_MINUTES", "15"))
if (is.na(IDLE_MINUTES) || IDLE_MINUTES < 1) IDLE_MINUTES <- 15

# --- Dataset registry -------------------------------------------------------
# key -> display label (order = navbar toggle order). Pancreatic first.
DATASET_LABELS <- c(
  mouse_pancreatic = "Pancreatic — KPP (mouse)",
  mouse_colorectal = "Colorectal — C26 (mouse)",
  human_pdac       = "PDAC (human)"
)

# Choose the contrast "tumour" arm: prefer the cachectic group (human), else the
# single non-Control condition (mouse KPP / C26). Must match data-raw/02 + 04.
pick_tumor <- function(conds) {
  nc <- setdiff(conds, "Control")
  if ("Cachectic" %in% nc) "Cachectic" else if (length(nc)) sort(nc)[1] else "Tumor"
}

.load_dataset <- function(k) {
  d <- file.path(app_data_dir, k)
  meta <- readRDS(file.path(d, "umap_metadata.rds"))
  tumor <- pick_tumor(unique(meta$condition))
  list(
    key          = k,
    label        = unname(DATASET_LABELS[[k]]),
    meta         = meta,
    gene_list    = readRDS(file.path(d, "gene_list.rds")),
    palette      = readRDS(file.path(d, "cell_type_palette.rds")),
    summary      = readRDS(file.path(d, "summary_stats.rds")),
    de_results   = readRDS(file.path(d, "de_results.rds")),
    marker_genes = readRDS(file.path(d, "marker_genes.rds")),
    # Pathway artifacts are nullable — the app still boots on a partial build
    # (before 04_pathways.R has run) and the Pathway tab shows an empty state.
    gsea         = if (file.exists(file.path(d, "pathways.rds"))) readRDS(file.path(d, "pathways.rds")) else NULL,
    progeny      = if (file.exists(file.path(d, "progeny.rds")))  readRDS(file.path(d, "progeny.rds"))  else NULL,
    cell_types   = sort(unique(meta$cell_type)),
    conditions   = sort(unique(meta$condition)),
    tumor_label  = tumor
  )
}

# Only register datasets whose artifacts actually exist (robust to partial builds).
AVAILABLE <- names(DATASET_LABELS)[
  vapply(names(DATASET_LABELS), function(k)
    file.exists(file.path(app_data_dir, k, "umap_metadata.rds")) &&
    file.exists(file.path(app_data_dir, k, "expression.h5")),
    logical(1))
]
DS <- lapply(AVAILABLE, .load_dataset)
names(DS) <- AVAILABLE
DEFAULT_DATASET <- AVAILABLE[1]

# --- Source helpers, open expression stores, then modules -------------------
for (f in list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)) source(f)
init_expr_stores(AVAILABLE)
for (f in list.files("modules", pattern = "\\.[Rr]$", full.names = TRUE)) source(f)
