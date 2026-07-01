#!/usr/bin/env Rscript
# 04_pathways.R — precompute pathway enrichment for the Pathway Enrichment tab.
#
# Per dataset, writes two tidy artifacts the app loads (no fgsea/progeny at runtime):
#   app/data/<ds>/pathways.rds   GSEA (fgsea) per cell_type x collection
#   app/data/<ds>/progeny.rds    PROGENy 14-pathway activity contrast per cell_type
#
# GSEA: rank each cell type's Tumor-vs-Control DE by the scanpy Wilcoxon z
# (`score` column — finite under p underflow, ~no ties), run fgsea (multilevel)
# against MSigDB Hallmark / Reactome / GO:BP, BH within (collection x cell_type).
# Sign convention everywhere: NES/contrast > 0 = UP in tumor/cachexia.
#
# PROGENy: per-cell activity (progeny, top=100), z-scored per pathway across the
# dataset; contrast = mean(z|tumor) - mean(z|control) per cell type (+ Cohen's d);
# a pseudobulk (per-sample) t-test supplies an exploratory p-value where >=2
# samples/condition exist (Neg/Pos sorts are NOT independent replicates, so these
# p-values are exploratory — the effect size is the headline).
#
# Run (base R 4.5.3, after 03_build_app_data.R):  Rscript data-raw/04_pathways.R

suppressPackageStartupMessages({
  library(fgsea); library(progeny); library(data.table)
  library(rhdf5); library(msigdbr)
})

.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
REPO <- if (length(.file_arg) == 1L) {
  dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
} else normalizePath(getwd())
APP_DATA <- file.path(REPO, "app", "data")
MSIG_MM  <- "/Users/calebgammon/rnaseq/c26_timecourse/published_style/data/processed/msigdbr_mm.rds"

DATASETS <- c(mouse_pancreatic = "Mouse", mouse_colorectal = "Mouse", human_pdac = "Human")
COLLECTIONS <- c("H", "C2_CP_REACTOME", "C5_GO_BP")

MINSIZE <- 15; MAXSIZE <- 500; SEED <- 0x5EED
PROGENY_TOP <- 100; MIN_GROUP_CELLS <- 10

# Contrast "tumour" arm: prefer cachectic (human), else the single non-Control
# condition (mouse). Must match data-raw/02 + app/global.R.
pick_tumor <- function(conds) {
  nc <- setdiff(conds, "Control")
  if ("Cachectic" %in% nc) "Cachectic" else if (length(nc)) sort(nc)[1] else NA_character_
}

# ---- gene-set sourcing -----------------------------------------------------
gene_sets_mouse <- function() {
  x <- readRDS(MSIG_MM)              # list(H, C2_CP_REACTOME, C5_GO_BP) -> named gene lists
  x[COLLECTIONS]
}
.msig <- function(coll, sub) {
  f <- function() if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = coll)
                  else msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub)
  g <- function() if (is.null(sub)) msigdbr(species = "Homo sapiens", category = coll)
                  else msigdbr(species = "Homo sapiens", category = coll, subcategory = sub)
  df <- tryCatch(f(), error = function(e) g())
  split(df$gene_symbol, df$gs_name)
}
gene_sets_human <- function() list(
  H              = .msig("H",  NULL),
  C2_CP_REACTOME = .msig("C2", "CP:REACTOME"),
  C5_GO_BP       = .msig("C5", "GO:BP")
)

pretty_label <- function(x) {
  x <- sub("^HALLMARK_", "", x); x <- sub("^GOBP_", "", x); x <- sub("^REACTOME_", "", x)
  tools::toTitleCase(tolower(gsub("_", " ", x)))
}

# ---- GSEA ------------------------------------------------------------------
rank_stat <- function(d) {
  d <- d[is.finite(d$score) & is.finite(d$avg_log2FC), ]
  d <- d[!duplicated(d$gene), ]
  stat <- d$score + sign(d$score) * abs(d$avg_log2FC) * 1e-9   # deterministic tie-break
  names(stat) <- d$gene
  sort(stat, decreasing = TRUE)
}

run_gsea <- function(de, sets, tumor_label) {
  out <- list(); skipped <- c()
  for (ct in names(de)) {
    stat <- rank_stat(de[[ct]])
    if (length(stat) < 200 || sum(de[[ct]]$p_val < 0.05, na.rm = TRUE) < 50) {
      skipped[ct] <- sprintf("too few informative genes (n=%d)", length(stat)); next
    }
    for (coll in COLLECTIONS) {
      set.seed(SEED)
      res <- tryCatch(
        fgsea::fgsea(pathways = sets[[coll]], stats = stat, minSize = MINSIZE,
                     maxSize = MAXSIZE, eps = 0, scoreType = "std",
                     nPermSimple = 1000, BPPARAM = BiocParallel::SerialParam()),
        error = function(e) NULL)
      if (is.null(res) || !nrow(res)) next
      le <- res$leadingEdge
      out[[length(out) + 1L]] <- data.table(
        collection = coll, cell_type = ct, pathway = res$pathway,
        pathway_label = pretty_label(res$pathway), size = res$size,
        ES = res$ES, NES = res$NES, pval = res$pval, padj = res$padj,
        direction = ifelse(res$NES >= 0, paste("Up in", tumor_label),
                           paste("Down in", tumor_label)),
        leadingEdge = vapply(le, function(g) paste(utils::head(g, 25), collapse = ","), character(1)),
        leadingEdge_n = lengths(le))
    }
  }
  gsea <- if (length(out)) rbindlist(out) else data.table()
  if (nrow(gsea)) {
    gsea[, collection := factor(collection, levels = COLLECTIONS)]
    gsea[, cell_type := factor(cell_type)]
    gsea[, direction := factor(direction)]
  }
  list(gsea = as.data.frame(gsea), skipped = skipped)
}

# ---- PROGENy ---------------------------------------------------------------
read_matrix <- function(h5, genes_present) {
  allg <- as.character(rhdf5::h5read(h5, "genes"))
  idx  <- match(genes_present, allg)
  M <- rhdf5::h5read(h5, "expr", index = list(NULL, idx))   # cells x genes (rhdf5 transposes)
  colnames(M) <- genes_present
  t(M)                                                       # genes x cells for progeny
}

run_progeny <- function(ds, organism, meta, gene_list) {
  mod_genes <- rownames(progeny::getModel(organism, top = PROGENY_TOP))
  present <- intersect(mod_genes, gene_list)
  h5 <- file.path(APP_DATA, ds, "expression.h5")
  mat <- read_matrix(h5, present)                            # genes x cells
  act <- progeny::progeny(mat, scale = FALSE, organism = organism,
                          top = PROGENY_TOP, perm = 1)        # cells x 14
  act <- scale(act)                                          # z per pathway across all cells
  act[is.nan(act)] <- 0

  conds <- meta$condition; tumor <- pick_tumor(unique(conds))
  samp <- meta$sample_id
  paths <- colnames(act); rows <- list()
  for (ct in sort(unique(meta$cell_type))) {
    cix <- meta$cell_type == ct
    t_i <- cix & conds == tumor; c_i <- cix & conds == "Control"
    if (sum(t_i) < MIN_GROUP_CELLS || sum(c_i) < MIN_GROUP_CELLS) next
    # per-sample pseudobulk means for the exploratory t-test
    samp_ct <- samp[cix]; act_ct <- act[cix, , drop = FALSE]; cond_ct <- conds[cix]
    for (p in paths) {
      a <- act[, p]
      contrast <- mean(a[t_i]) - mean(a[c_i])
      pooled_sd <- stats::sd(a[cix]); cohens_d <- if (pooled_sd > 0) contrast / pooled_sd else NA_real_
      # pseudobulk: mean activity per sample, split by condition
      pb <- tapply(act_ct[, p], samp_ct, mean)
      pb_cond <- tapply(cond_ct, samp_ct, function(x) x[1])[names(pb)]
      tv <- pb[pb_cond == tumor]; cv <- pb[pb_cond == "Control"]
      low_n <- length(tv) < 2 || length(cv) < 2
      tt <- if (low_n) NULL else tryCatch(stats::t.test(tv, cv), error = function(e) NULL)
      rows[[length(rows) + 1L]] <- data.table(
        cell_type = ct, pathway = p, act_tumor = mean(a[t_i]), act_control = mean(a[c_i]),
        contrast = contrast, cohens_d = cohens_d,
        t_stat = if (is.null(tt)) NA_real_ else unname(tt$statistic),
        pval = if (is.null(tt)) NA_real_ else tt$p.value,
        n_tumor = sum(t_i), n_control = sum(c_i),
        n_samp_tumor = length(tv), n_samp_control = length(cv), low_n = low_n)
    }
  }
  pr <- rbindlist(rows)
  pr[, padj := p.adjust(pval, "BH")]
  pr[, cell_type := factor(cell_type)]; pr[, pathway := factor(pathway)]
  attr(pr, "tumor_label") <- tumor
  as.data.frame(pr)
}

# ---- biological sanity checks ---------------------------------------------
sanity <- function(ds, gsea, progeny) {
  cat(sprintf("[04] --- sanity: %s ---\n", ds))
  g <- as.data.table(gsea)
  chk <- function(path, ct, expect_up = TRUE) {
    r <- g[pathway == path & cell_type == ct]
    if (!nrow(r)) { cat(sprintf("   ?  %-38s %-20s (not tested)\n", path, ct)); return(NA) }
    ok <- if (expect_up) r$NES > 0 else r$NES < 0
    cat(sprintf("   %s %-38s %-18s NES=%+.2f padj=%.1e\n", ifelse(ok, "PASS", "FAIL"),
                path, ct, r$NES, r$padj)); ok
  }
  mac <- "Macrophage/Monocyte"
  chk("HALLMARK_TNFA_SIGNALING_VIA_NFKB", mac, TRUE)
  chk("HALLMARK_INFLAMMATORY_RESPONSE", mac, TRUE)
  chk("HALLMARK_INTERFERON_GAMMA_RESPONSE", mac, TRUE)
  pr <- as.data.table(progeny)
  for (p in c("NFkB", "TNFa")) {
    r <- pr[pathway == p & cell_type == mac]
    if (nrow(r)) cat(sprintf("   %s PROGENy %-6s %-18s contrast=%+.2f\n",
                             ifelse(r$contrast > 0, "PASS", "FAIL"), p, mac, r$contrast))
  }
}

# ---- driver ----------------------------------------------------------------
# Optional: Rscript data-raw/04_pathways.R --dataset human_pdac  (re-run one)
.only <- grep("^--dataset=", commandArgs(TRUE), value = TRUE)
.only <- if (length(.only)) sub("^--dataset=", "", .only[1]) else
         { a <- commandArgs(TRUE); i <- match("--dataset", a); if (!is.na(i)) a[i + 1] else NULL }
TARGETS <- if (!is.null(.only)) intersect(.only, names(DATASETS)) else names(DATASETS)

for (ds in TARGETS) {
  organism <- DATASETS[[ds]]
  cat(sprintf("\n[04] ==== %s (%s) ====\n", ds, organism))
  d <- file.path(APP_DATA, ds)
  if (!file.exists(file.path(d, "de_results.rds"))) { cat("   skip (no de_results)\n"); next }
  meta <- readRDS(file.path(d, "umap_metadata.rds"))
  de   <- readRDS(file.path(d, "de_results.rds"))
  gene_list <- readRDS(file.path(d, "gene_list.rds"))
  tumor <- pick_tumor(unique(meta$condition))
  sets <- if (organism == "Mouse") gene_sets_mouse() else gene_sets_human()
  cat(sprintf("[04]   gene sets: H=%d Reactome=%d GO:BP=%d | Hallmark coverage=%.0f%%\n",
              length(sets$H), length(sets$C2_CP_REACTOME), length(sets$C5_GO_BP),
              100 * mean(unique(unlist(sets$H)) %in% gene_list)))

  cat("[04]   running GSEA (fgsea multilevel)...\n")
  gres <- run_gsea(de, sets, tumor)
  attr(gres$gsea, "params") <- list(ranking = "scanpy Wilcoxon z (score)",
    minSize = MINSIZE, maxSize = MAXSIZE, eps = 0, scoreType = "std", seed = SEED,
    fgsea = as.character(utils::packageVersion("fgsea")),
    msigdbr = as.character(utils::packageVersion("msigdbr")), skipped = gres$skipped)
  saveRDS(gres$gsea, file.path(d, "pathways.rds"), compress = "xz")
  cat(sprintf("[04]   pathways.rds: %d rows (%d cell types tested)\n",
              nrow(gres$gsea), length(unique(gres$gsea$cell_type))))

  cat("[04]   running PROGENy...\n")
  pr <- run_progeny(ds, organism, meta, gene_list)
  attr(pr, "params") <- list(organism = organism, top = PROGENY_TOP,
    scaling = "dataset z-score per pathway", inference = "pseudobulk t-test",
    progeny = as.character(utils::packageVersion("progeny")))
  saveRDS(pr, file.path(d, "progeny.rds"), compress = "xz")
  cat(sprintf("[04]   progeny.rds: %d rows (%d cell types)\n", nrow(pr), length(unique(pr$cell_type))))

  sanity(ds, gres$gsea, pr)
}
cat("\n[04] Done.\n")
