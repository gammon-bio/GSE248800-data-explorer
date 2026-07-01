"""
02_export_from_h5ad.py — build the Shiny app's heavy artifacts from the three
annotated AnnData objects written by 01_process.py.

For each dataset (mouse_pancreatic, mouse_colorectal, human_pdac) this writes
into app/data/<dataset>/:

  expression.h5              HDF5 (genes x cells, log-norm, (1,n_cells) chunks,
                             gzip-4) — the EXACT format app/R/expr_store.R reads.
  de_<celltype>.csv          Tumor-vs-Control Wilcoxon DE within each cell type,
                             columns: gene, avg_log2FC, p_val, p_val_adj,
                             pct.1 (tumor), pct.2 (control), cell_type.
  marker_genes.csv           one-vs-rest cell-type markers (gene, cluster,
                             avg_log2FC, pct.1, pct.2) — top genes per cell type.
  _scaffold/umap_metadata.csv, genes.txt   consumed by 03_build_app_data.R.

Invariants (mirrors kpc_tme): expr row order == genes, column order == cells ==
umap_metadata.csv row order; .X is LOG-NORM.

Run in the scanpy env:
  /opt/homebrew/Caskroom/miniforge/base/envs/zhang_atlas/bin/python data-raw/02_export_from_h5ad.py
"""

import os
import re

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR = os.path.join(REPO, "data-raw", "data")
APP_DATA = os.path.join(REPO, "app", "data")

DATASETS = ["mouse_pancreatic", "mouse_colorectal", "human_pdac"]

GENE_BLOCK = 2000       # gene rows per densify block when writing expression.h5
MIN_CELLS_DE = 10       # per-side minimum for a Tumor-vs-Control cell-type test
TOP_MARKERS = 50        # markers kept per cell type in marker_genes.csv

# Human PDAC has three groups; the biologically meaningful cachexia contrast is
# CACHECTIC vs Control (weight-stable PDAC dilutes the signal). We therefore use
# the fine condition_detail as the app-facing `condition`, relabelled cleanly,
# and pick the cachectic arm as "tumour". Weight-stable cells stay in the atlas
# (UMAP / composition) but are excluded from the Cachectic-vs-Control contrast
# (rank_genes_groups with reference="Control", groups=[Cachectic] ignores them).
CLEAN_COND = {"PDAC_cachectic": "Cachectic", "PDAC_weight_stable": "Weight-stable"}


def clean_cond(x):
    return CLEAN_COND.get(str(x), str(x))


def pick_tumor(conds):
    """Choose the contrast 'tumour' arm: prefer the cachectic group, else the
    single non-Control condition (mouse KPP / C26)."""
    nc = [c for c in conds if c != "Control"]
    if "Cachectic" in nc:
        return "Cachectic"
    return sorted(nc)[0] if nc else None


def _san(name):
    """Filesystem-safe cell-type token for de_<...>.csv filenames."""
    return re.sub(r"[^0-9A-Za-z]+", "_", str(name)).strip("_")


def _write_expression_h5(adata, path):
    """genes x cells, float32, log-norm, per-gene chunks, gzip-4 (kpc_tme format).
    Idempotent: the gene/cell order is deterministic from the annotated h5ad, so
    if expression.h5 already exists we return that order and skip the (expensive)
    rewrite — lets us cheaply re-export just the DE/marker CSVs."""
    n_genes, n_cells = adata.n_vars, adata.n_obs
    genes = adata.var_names.to_numpy().astype(str)
    cells = adata.obs_names.to_numpy().astype(str)
    if os.path.exists(path):
        print(f"[02]   expression.h5 exists ({os.path.getsize(path)/1e6:.0f} MB) — skipping rewrite")
        return genes, cells
    X = adata.X.tocsc() if sparse.issparse(adata.X) else np.asarray(adata.X)
    str_dt = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as h5:
        expr = h5.create_dataset(
            "expr", shape=(n_genes, n_cells), dtype="float32",
            chunks=(1, n_cells), compression="gzip", compression_opts=4,
        )
        for start in range(0, n_genes, GENE_BLOCK):
            stop = min(start + GENE_BLOCK, n_genes)
            block = X[:, start:stop]
            if sparse.issparse(block):
                block = block.toarray()
            expr[start:stop, :] = np.asarray(block, dtype="float32").T
        h5.create_dataset("genes", data=np.array(genes, dtype=object), dtype=str_dt)
        h5.create_dataset("cells", data=np.array(cells, dtype=object), dtype=str_dt)
    return genes, cells


def _pct_expressing(adata, mask):
    """% of cells in `mask` with X>0, per gene (var order)."""
    sub = adata.X[mask]
    nz = np.asarray((sub > 0).sum(axis=0)).ravel().astype(float)
    denom = max(int(mask.sum()), 1)
    return nz / denom * 100.0


def _de_tumor_vs_control(adata, genes, out_dir):
    """Per cell type: Wilcoxon tumor vs Control. Returns list of written cell types."""
    conds = adata.obs["condition"].astype(str)
    tumor = pick_tumor(conds.unique())   # Cachectic (human) / KPP / C26
    if tumor is None:
        print("[02]   WARN no non-Control condition; skipping DE")
        return []
    written = []
    for ct in sorted(adata.obs["cell_type"].astype(str).unique()):
        sub = adata[adata.obs["cell_type"].astype(str) == ct].copy()
        sc_ = sub.obs["condition"].astype(str)
        n_t, n_c = int((sc_ == tumor).sum()), int((sc_ == "Control").sum())
        if n_t < MIN_CELLS_DE or n_c < MIN_CELLS_DE:
            print(f"[02]     skip {ct}: tumor={n_t} control={n_c} (< {MIN_CELLS_DE})")
            continue
        sc.tl.rank_genes_groups(sub, "condition", groups=[tumor],
                                reference="Control", method="wilcoxon")
        df = sc.get.rank_genes_groups_df(sub, group=tumor).set_index("names")
        pct_t = pd.Series(_pct_expressing(sub, (sc_ == tumor).to_numpy()), index=sub.var_names)
        pct_c = pd.Series(_pct_expressing(sub, (sc_ == "Control").to_numpy()), index=sub.var_names)
        out = pd.DataFrame({
            "gene": genes,
            "avg_log2FC": df["logfoldchanges"].reindex(genes).to_numpy(),
            # Wilcoxon z-statistic (signed, continuous) — the GSEA ranking metric.
            # Unlike -log10(p) it does not saturate at the p=0 underflow ceiling,
            # so it avoids the massive top-of-list ties a p-based rank produces.
            "score": df["scores"].reindex(genes).to_numpy(),
            "p_val": df["pvals"].reindex(genes).to_numpy(),
            "p_val_adj": df["pvals_adj"].reindex(genes).to_numpy(),
            "pct.1": pct_t.reindex(genes).to_numpy(),   # tumor
            "pct.2": pct_c.reindex(genes).to_numpy(),   # control
            "cell_type": ct,
        })
        out = out[np.isfinite(out["avg_log2FC"])]
        out.to_csv(os.path.join(out_dir, f"de_{_san(ct)}.csv"), index=False)
        written.append(ct)
        print(f"[02]     de {ct}: tumor={n_t} vs control={n_c} ({tumor} vs Control)")
    return written


def _markers_one_vs_rest(adata, genes, out_path):
    """One-vs-rest Wilcoxon per cell type; keep TOP_MARKERS positive-LFC genes each."""
    sc.tl.rank_genes_groups(adata, "cell_type", method="wilcoxon")
    rows = []
    for ct in adata.obs["cell_type"].astype(str).unique():
        df = sc.get.rank_genes_groups_df(adata, group=ct)
        df = df[df["logfoldchanges"] > 0].sort_values("logfoldchanges", ascending=False).head(TOP_MARKERS)
        rows.append(pd.DataFrame({
            "gene": df["names"].to_numpy(),
            "cluster": ct,
            "avg_log2FC": df["logfoldchanges"].to_numpy(),
            "p_val_adj": df["pvals_adj"].to_numpy(),
        }))
    markers = pd.concat(rows, ignore_index=True)
    markers.to_csv(out_path, index=False)
    print(f"[02]   marker_genes.csv: {len(markers)} rows across "
          f"{adata.obs['cell_type'].nunique()} cell types")


def export(dataset):
    src = os.path.join(SRC_DIR, f"{dataset}_annotated.h5ad")
    if not os.path.exists(src):
        print(f"[02]   SKIP {dataset}: {src} not found")
        return
    print(f"\n[02] ==== {dataset} ====")
    adata = ad.read_h5ad(src)
    print(f"[02]   {adata.n_obs} cells x {adata.n_vars} genes, "
          f"{adata.obs['cell_type'].nunique()} cell types")

    # App-facing `condition` = the fine group, cleanly relabelled. For mouse this
    # is unchanged (KPP/C26/Control); for human it becomes the 3-level
    # Control / Cachectic / Weight-stable so the contrast can be Cachectic-vs-Control.
    adata.obs["condition"] = adata.obs["condition_detail"].astype(str).map(clean_cond)

    out_dir = os.path.join(APP_DATA, dataset)
    scaffold = os.path.join(out_dir, "_scaffold")
    os.makedirs(scaffold, exist_ok=True)

    # Remove stale de_<ct>.csv so a changed contrast (e.g. switching human to
    # Cachectic-vs-Control, which skips T/NK) can't leave an old table behind.
    import glob as _glob
    for f in _glob.glob(os.path.join(out_dir, "de_*.csv")):
        os.remove(f)

    genes, cells = _write_expression_h5(adata, os.path.join(out_dir, "expression.h5"))
    print(f"[02]   expression.h5: ({len(genes)}, {len(cells)}) float32 gzip-4")

    umap = adata.obsm["X_umap"]
    meta = pd.DataFrame({
        "cell": cells,
        "UMAP1": np.asarray(umap[:, 0], dtype=float),
        "UMAP2": np.asarray(umap[:, 1], dtype=float),
        "cell_type": adata.obs["cell_type"].astype(str).to_numpy(),
        "condition": adata.obs["condition"].astype(str).to_numpy(),
        "condition_detail": adata.obs["condition_detail"].astype(str).to_numpy(),
        "sample_id": adata.obs["sample_id"].astype(str).to_numpy(),
        "dataset": dataset,
    })
    meta.to_csv(os.path.join(scaffold, "umap_metadata.csv"), index=False)
    with open(os.path.join(scaffold, "genes.txt"), "w") as f:
        f.write("\n".join(genes.tolist()) + "\n")

    _de_tumor_vs_control(adata, genes, out_dir)
    _markers_one_vs_rest(adata, genes, os.path.join(out_dir, "marker_genes.csv"))
    print(f"[02]   done {dataset}")


if __name__ == "__main__":
    for ds in DATASETS:
        export(ds)
    print("\n[02] Done. Next: Rscript data-raw/03_build_app_data.R")
