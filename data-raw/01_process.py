"""
01_process.py — scRNA-seq processing for the GSE248800 non-myofiber explorer.

GSE248800 ("Muscle Inflammation is Regulated by NF-kB ... in Cancer Cachexia")
is whole-cell scRNA-seq of skeletal-muscle mononuclear cells. This script turns
the raw per-sample 10x H5 files (downloaded by 00_download.sh) into three
independently processed, annotated AnnData objects — one per dataset the Shiny
app toggles between:

  mouse_pancreatic  (KPP model,  GSM7920260-263)
  mouse_colorectal  (C26 model,  GSM7920264-270)
  human_pdac        (human PDAC, GSM7920254-259)

For each dataset: read + concatenate samples -> QC filter -> normalize/log1p ->
HVG/PCA/neighbors/Leiden/UMAP -> marker-score annotation of NON-MYOFIBER
lineages -> DROP myofiber/contaminant clusters -> write
data-raw/data/<dataset>_annotated.h5ad.

The NF-kB reporter Neg/Pos fractions are POOLED per condition (they are just two
sorted fractions of the same tissue). The primary contrast the app exposes is
Tumor-vs-Control, carried in obs['condition'] (2-level) with obs['condition_detail']
keeping the human weight-stable/cachectic split.

Invariants downstream code depends on:
  - adata.X is LOG-NORMALIZED (normalize_total 1e4 + log1p); raw counts live in
    layers['counts'].
  - obs has: cell_type, condition, condition_detail, sample_id, gsm, leiden.
  - No cell_type == "Myofiber" survives (those clusters are dropped).

Run in the scanpy env:
  /opt/homebrew/Caskroom/miniforge/base/envs/zhang_atlas/bin/python data-raw/01_process.py
Optional: --dataset mouse_pancreatic  (process just one)
"""

import argparse
import glob
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 1

# ----------------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------------
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(REPO, "data-raw", "data", "raw")
OUT_DIR = os.path.join(REPO, "data-raw", "data")
os.makedirs(OUT_DIR, exist_ok=True)

# GSM -> (dataset, tumor-vs-control condition, fine condition_detail).
# Titles come from the GEO record; Neg/Pos (NF-kB reporter fractions) are pooled.
SAMPLE_MAP = {
    # human PDAC muscle
    "GSM7920254": ("human_pdac", "PDAC", "PDAC_weight_stable"),   # Neg_WS_PDAC
    "GSM7920255": ("human_pdac", "PDAC", "PDAC_weight_stable"),   # Pos_WS_PDAC
    "GSM7920256": ("human_pdac", "Control", "Control"),           # Neg_Control
    "GSM7920257": ("human_pdac", "Control", "Control"),           # Pos_Control
    "GSM7920258": ("human_pdac", "PDAC", "PDAC_cachectic"),       # Neg_C_PDAC
    "GSM7920259": ("human_pdac", "PDAC", "PDAC_cachectic"),       # Pos_C_PDAC
    # mouse KPP (pancreatic)
    "GSM7920260": ("mouse_pancreatic", "Control", "Control"),     # Ctrl_Neg
    "GSM7920261": ("mouse_pancreatic", "KPP", "KPP"),             # Kpp_Neg
    "GSM7920262": ("mouse_pancreatic", "Control", "Control"),     # Ctrl_Pos
    "GSM7920263": ("mouse_pancreatic", "KPP", "KPP"),             # Kpp_Pos
    # mouse C26 (colorectal)
    "GSM7920264": ("mouse_colorectal", "Control", "Control"),     # negControl_s1
    "GSM7920265": ("mouse_colorectal", "C26", "C26"),            # C26neg_s1
    "GSM7920266": ("mouse_colorectal", "Control", "Control"),     # posControl_s1
    "GSM7920267": ("mouse_colorectal", "C26", "C26"),            # C26pos_s1
    "GSM7920268": ("mouse_colorectal", "Control", "Control"),     # posControl_s2
    "GSM7920269": ("mouse_colorectal", "C26", "C26"),            # C26pos_s2
    "GSM7920270": ("mouse_colorectal", "C26", "C26"),            # C26neg_s2
}

DATASET_ORGANISM = {
    "mouse_pancreatic": "mouse",
    "mouse_colorectal": "mouse",
    "human_pdac": "human",
}

# Non-myofiber lineage markers (mouse symbols; upper-cased for human).
# Cluster -> lineage assignment is by argmax of per-cluster mean marker score.
MARKERS = {
    "Macrophage/Monocyte": ["Ptprc", "Adgre1", "Cd68", "Lyz2", "Csf1r", "Itgam", "Fcgr1"],
    "Dendritic cell": ["Flt3", "Xcr1", "Cd209a", "Clec9a", "Batf3"],
    "T/NK cell": ["Cd3e", "Cd3d", "Cd8a", "Cd4", "Nkg7", "Klrb1c", "Gzmb"],
    "B cell": ["Cd19", "Ms4a1", "Cd79a", "Cd79b"],
    "Neutrophil": ["S100a8", "S100a9", "Ly6g", "Retnlg"],
    "FAP/Fibroblast": ["Pdgfra", "Dcn", "Col1a1", "Col3a1", "Lum", "Gsn"],
    "Endothelial": ["Pecam1", "Cdh5", "Cldn5", "Flt1", "Egfl7"],
    "Pericyte/Mural": ["Rgs5", "Pdgfrb", "Notch3", "Kcnj8", "Acta2"],
    "Satellite cell": ["Pax7", "Myf5", "Chodl"],
    "Adipocyte": ["Adipoq", "Plin1", "Lep", "Fabp4"],
    "Schwann/Glia": ["Plp1", "Mpz", "Sox10"],
    "Tenocyte": ["Scx", "Tnmd", "Mkx"],
}
# Myofiber / mature-muscle contamination — clusters scoring highest here are dropped.
MYOFIBER_MARKERS = ["Acta1", "Myh1", "Myh2", "Myh4", "Ttn", "Ckm", "Tnnt3", "Des", "Mb", "Tnnc2"]

# QC thresholds. Only RAW (unfiltered) 10x matrices were released, so every
# sample carries ~1M barcodes that are mostly empty droplets. min_genes acts as
# the empty-droplet knee and is applied PER SAMPLE before concatenation to bound
# memory; the remaining thresholds run on the merged object.
MIN_GENES = 200        # empty-droplet floor (per-sample, pre-concat)
MIN_COUNTS = 500       # total-UMI floor (drops lowest-quality barcodes)
MIN_CELLS = 3
MAX_PCT_MT = 20.0
MAX_GENES = 8000       # crude high-complexity/doublet ceiling
LEIDEN_RES = 1.0
N_HVG = 2000
N_PCS = 30


def _symbols_for(organism, symbols):
    """Mouse markers as-is; human -> upper-case (HGNC convention)."""
    return [s.upper() for s in symbols] if organism == "human" else list(symbols)


def _read_sample(h5_path):
    """Read one sample as AnnData, tolerating 10x H5 or generic h5ad."""
    try:
        a = sc.read_10x_h5(h5_path)
    except Exception:
        a = sc.read_h5ad(h5_path)
    a.var_names_make_unique()
    # Keep only Gene Expression features if a multi-feature 10x matrix.
    if "feature_types" in a.var.columns:
        ge = a.var["feature_types"].astype(str).str.contains("Gene Expression")
        if ge.any():
            a = a[:, ge.values].copy()
    # Empty-droplet knee: drop the ~1M near-empty barcodes now (pre-concat) so
    # only real cells are carried into the merge.
    sc.pp.filter_cells(a, min_genes=MIN_GENES)
    return a


def _load_dataset(dataset):
    """Concatenate every GSM belonging to `dataset` with obs annotations."""
    gsms = [g for g, (d, _, _) in SAMPLE_MAP.items() if d == dataset]
    parts = []
    for gsm in sorted(gsms):
        matches = glob.glob(os.path.join(RAW_DIR, dataset, f"{gsm}*")) + \
                  glob.glob(os.path.join(RAW_DIR, f"{gsm}*"))
        matches = [m for m in matches if m.endswith((".h5", ".h5ad"))]
        if not matches:
            print(f"[01]   WARN no H5 found for {gsm} in {dataset}")
            continue
        h5 = sorted(matches)[0]
        a = _read_sample(h5)
        _, cond, detail = SAMPLE_MAP[gsm]
        a.obs["gsm"] = gsm
        a.obs["sample_id"] = f"{gsm}"
        a.obs["condition"] = cond
        a.obs["condition_detail"] = detail
        a.obs_names = [f"{gsm}_{bc}" for bc in a.obs_names]
        print(f"[01]   {gsm}: {a.n_obs} cells x {a.n_vars} genes  ({os.path.basename(h5)})")
        parts.append(a)
    if not parts:
        raise SystemExit(f"No samples read for {dataset}")
    adata = ad.concat(parts, join="inner", index_unique=None)
    adata.obs_names_make_unique()
    return adata


def _qc_filter(adata, organism):
    mt_prefix = "MT-" if organism == "human" else "mt-"
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    sc.pp.filter_cells(adata, min_counts=MIN_COUNTS)
    adata = adata[adata.obs["n_genes_by_counts"] < MAX_GENES].copy()
    adata = adata[adata.obs["pct_counts_mt"] < MAX_PCT_MT].copy()
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    print(f"[01]   QC: {n0} -> {adata.n_obs} cells, {adata.n_vars} genes "
          f"(min_genes={MIN_GENES}, min_counts={MIN_COUNTS}, pct_mt<{MAX_PCT_MT}, n_genes<{MAX_GENES})")
    return adata


def _cluster(adata):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)  # adata.X is now LOG-NORM — this is what we export.

    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat")
    # PCA on a scaled HVG COPY so adata.X stays log-norm for the expression store.
    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=N_PCS, svd_solver="arpack")
    adata.obsm["X_pca"] = hvg.obsm["X_pca"]

    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=N_PCS, use_rep="X_pca")
    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph",
                     n_iterations=2, directed=False)
    except TypeError:  # older scanpy signature
        sc.tl.leiden(adata, resolution=LEIDEN_RES)
    sc.tl.umap(adata)
    print(f"[01]   clustered: {adata.obs['leiden'].nunique()} Leiden clusters")
    return adata


def _annotate(adata, organism):
    """Score each lineage, assign every Leiden cluster its argmax lineage.
    Clusters whose dominant score is myofiber are labeled 'Myofiber' for removal."""
    score_cols = {}
    for lineage, syms in {**MARKERS, "Myofiber": MYOFIBER_MARKERS}.items():
        genes = [g for g in _symbols_for(organism, syms) if g in adata.var_names]
        col = f"score_{lineage}"
        if genes:
            sc.tl.score_genes(adata, genes, score_name=col)
        else:
            adata.obs[col] = 0.0
        score_cols[lineage] = col

    # Per-cluster mean score -> argmax lineage.
    lineages = list(score_cols)
    means = adata.obs.groupby("leiden", observed=True)[list(score_cols.values())].mean()
    cluster_label = {}
    for cl, row in means.iterrows():
        best = lineages[int(np.argmax(row.values))]
        cluster_label[cl] = best
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_label).astype(str)

    # Report + drop myofiber clusters.
    labelled = adata.obs.groupby("cell_type", observed=True).size().sort_values(ascending=False)
    print("[01]   cluster labels:")
    for ct, n in labelled.items():
        print(f"[01]      {ct}: {n}")
    keep = adata.obs["cell_type"] != "Myofiber"
    n_drop = int((~keep).sum())
    adata = adata[keep].copy()
    print(f"[01]   dropped {n_drop} myofiber-cluster cells -> {adata.n_obs} non-myofiber cells")
    # Drop the transient score columns before saving.
    adata.obs.drop(columns=list(score_cols.values()), inplace=True, errors="ignore")
    return adata


def process(dataset):
    organism = DATASET_ORGANISM[dataset]
    print(f"\n[01] ==== {dataset} ({organism}) ====")
    adata = _load_dataset(dataset)
    print(f"[01]   merged: {adata.n_obs} cells x {adata.n_vars} genes")
    adata = _qc_filter(adata, organism)
    adata = _cluster(adata)
    adata = _annotate(adata, organism)

    # Tidy obs for export.
    adata.obs["dataset"] = dataset
    adata.obs["organism"] = organism
    keep_cols = ["cell_type", "condition", "condition_detail", "sample_id",
                 "gsm", "dataset", "organism", "leiden",
                 "n_genes_by_counts", "total_counts", "pct_counts_mt"]
    adata.obs = adata.obs[[c for c in keep_cols if c in adata.obs.columns]]

    out = os.path.join(OUT_DIR, f"{dataset}_annotated.h5ad")
    adata.write_h5ad(out)
    print(f"[01]   wrote {out}  ({adata.n_obs} cells, {adata.n_vars} genes)")
    print(f"[01]   conditions: {dict(adata.obs['condition'].value_counts())}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=list(DATASET_ORGANISM), default=None,
                    help="process a single dataset (default: all three)")
    args = ap.parse_args()
    datasets = [args.dataset] if args.dataset else list(DATASET_ORGANISM)
    for ds in datasets:
        process(ds)
    print("\n[01] Done. Next: 02_export_from_h5ad.py")
