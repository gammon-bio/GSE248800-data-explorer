# GSE248800 — Muscle Non-myofiber scRNA-seq Explorer

Interactive Shiny app for the **non-myofiber (mononuclear) compartment** of
skeletal muscle in cancer cachexia, from **GSE248800** (“Muscle Inflammation is
Regulated by NF-κB from Multiple Cells to Control Distinct States of Wasting in
Cancer Cachexia”). Whole-cell scRNA-seq captures the stromal/immune populations
— macrophages, FAPs/fibroblasts, endothelial, pericytes, T/NK/B cells,
neutrophils, satellite cells and more.

A **dataset toggle** in the navbar switches between three independently processed
datasets:

| Dataset | Model | Non-myofiber cells | Tumor vs Control |
|---|---|---|---|
| Pancreatic — KPP (mouse) | KPP GEMM (GSM7920260–263) | 50,567 | KPP vs Control |
| Colorectal — C26 (mouse) | C26 (GSM7920264–270) | 32,042 | C26 vs Control |
| PDAC (human) | human PDAC (GSM7920254–259) | 13,173 | PDAC vs Control |

## Tabs
- **Gene Explorer** — UMAP coloured by a gene's log-norm expression + violin/dot across cell types.
- **Differential Expression** — Tumor-vs-Control volcano + table, per cell type.
- **Pathway Enrichment** — top up/down pathways by cell type in cachexia: GSEA (Hallmark / Reactome /
  GO:BP) NES heatmap + top-pathway bars + a pathway-across-cell-types view, plus a **PROGENy** 14-pathway
  activity heatmap (NF-κB / TNFa / JAK-STAT …). Red = up in cachexia.
- **Cell Composition** — UMAP by cell type / by condition, proportions bar, counts.
- **Marker Genes** — top one-vs-rest marker heatmap + custom-gene dot plot.
- **About** — methods, citation, per-dataset summary.

Each tab uses sub-tabs (one plot at a time) and scrolls; every plot downloads as PDF/CSV.

## Run

The app runs on **base R 4.5.3** (needs `shiny`, `bslib`, `shinyWidgets`, `DT`,
`plotly`, `ggplot2`, `dplyr`, `tidyr`, `viridis`, `scales`, `waiter`, `rhdf5` —
all in the default library here; no Seurat at runtime).

```bash
Rscript -e 'shiny::runApp("app", port = 8765)'
# then open http://127.0.0.1:8765
```

## Deploy (Docker)

The app bundles ~675 MB of HDF5 expression data, so it ships as a self-contained
container (data baked in, no upload limits). Runtime needs only the shiny stack +
`rhdf5`; `fgsea`/`progeny`/`msigdbr` are build-time only.

```bash
docker build -t bryce-scrnaseq .            # from repo root (~5–10 min first build)
docker run --rm -p 3838:3838 bryce-scrnaseq # then open http://localhost:3838
```

The same image runs anywhere a container does:

- **Google Cloud Run** — `gcloud run deploy bryce-scrnaseq --source . --port 3838 --memory 2Gi --cpu 2 --allow-unauthenticated` (or push the built image to Artifact Registry and `--image`). Set `--min-instances 0` to scale to zero.
- **A VPS** — `docker run -d --restart unless-stopped -p 80:3838 bryce-scrnaseq` behind nginx/Caddy for TLS.
- **shinyapps.io** is not recommended here (the 675 MB data strains its bundle/tier limits); use a container instead.

Give the container ~2 GB RAM: the HDF5 matrices stay on disk (read one gene at a
time), so RAM is dominated by the metadata/DE/pathway tables, not the expression.

## Rebuild the data (`data-raw/`)

The raw data is downloaded from GEO and processed with **scanpy** (conda env
`zhang_atlas`). Heavy artifacts (`expression.h5`, `*_annotated.h5ad`, the raw
tar) are git-ignored — rebuild them with:

```bash
bash   data-raw/00_download.sh                    # GSE248800_RAW.tar -> 17 raw 10x H5
ENV=/opt/homebrew/Caskroom/miniforge/base/envs/zhang_atlas/bin/python
$ENV data-raw/01_process.py            # QC -> Leiden -> annotate -> drop myofibers -> *_annotated.h5ad
$ENV data-raw/02_export_from_h5ad.py   # -> app/data/<dataset>/expression.h5 + DE + marker CSVs
Rscript data-raw/03_build_app_data.R   # -> the .rds artifacts the app loads
Rscript data-raw/04_pathways.R         # -> pathways.rds (fgsea) + progeny.rds per dataset
```

`04_pathways.R` needs `fgsea`, `progeny`, `msigdbr` (installed in the base R library). GSEA ranks each
cell type's Tumor-vs-Control genes by the scanpy Wilcoxon z (`score` column) against MSigDB Hallmark /
Reactome / GO:BP (BH within collection × cell type); PROGENy scores 14 pathway activities per cell,
reported as the cachexia − control contrast (effect size + exploratory pseudobulk p). The pathway
artifacts are nullable — the app boots without them and shows an empty-state.

## Method notes / caveats
- Only **raw (unfiltered)** 10x matrices were released; an empty-droplet knee
  filter (≥200 genes / ≥500 UMIs, <20% mito, per sample) approximates cell calling.
- NF-κB reporter **Neg/Pos** fractions are **pooled** per condition.
- Each dataset is processed **independently** (no cross-dataset integration).
- Cell-type labels are **marker-score defaults** for exploration and should be
  validated experimentally.
- DE is per-cell Wilcoxon (tumor vs control); with thousands of cells most genes
  clear p<0.05, so use the log2FC / −log10(p) sliders to focus on strong effects.
- **Pathway sign:** NES / PROGENy contrast > 0 = up in cachexia (red). The NF-κB /
  TNFa / inflammatory signature is robustly up in **mouse** (KPP, C26) myeloid cells.
- **Human contrast = cachectic vs control:** the human contrast targets true cachexia
  (weight-stable PDAC cells stay in the atlas but are excluded from DE/pathways; T/NK is
  dropped for too few cachectic cells). This recovers signal in well-powered cell types
  (e.g. FAP/Fibroblast), but human macrophages have only ~10 cachectic cells so remain
  underpowered — the mouse models carry the robust myeloid inflammation story.
- PROGENy p-values are pseudobulk/sample-level and exploratory (NF-κB reporter
  Neg/Pos fractions are not independent replicates); the effect size is the headline.

Architecture reuses the HDF5-backed engine from `kpc_tme` (lazy per-gene reads,
so the ~700 MB of expression never loads into memory) and the tab set from
`zhang_shiny`.
