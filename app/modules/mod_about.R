# mod_about.R — static description + methods + a dataset-aware summary table.

aboutUI <- function(id) {
  ns <- NS(id)
  layout_column_wrap(
    width = 1, fill = FALSE,
    card(
      card_header("About this app"),
      card_body(
        h4("GSE248800 — Muscle non-myofiber scRNA-seq explorer"),
        p("Interactive explorer of the ", strong("non-myofiber (mononuclear)"),
          " compartment of skeletal muscle in cancer cachexia, from ",
          strong("GSE248800"), " (“Muscle Inflammation is Regulated by NF-κB ",
          "from Multiple Cells to Control Distinct States of Wasting in Cancer ",
          "Cachexia”). Whole-cell scRNA-seq captures the stromal and immune ",
          "populations — macrophages, FAPs/fibroblasts, endothelial, pericytes, ",
          "T/NK/B cells and more — that drive muscle inflammation."),
        p("Use the ", strong("dataset toggle"), " in the top bar to switch between:"),
        tags$ul(
          tags$li(strong("Pancreatic — KPP (mouse):"), " KPP tumour-bearing vs control muscle."),
          tags$li(strong("Colorectal — C26 (mouse):"), " C26 tumour-bearing vs control muscle."),
          tags$li(strong("PDAC (human):"), " human PDAC muscle — the contrast is ",
                  strong("cachectic PDAC vs control"), " (weight-stable PDAC cells are ",
                  "shown in the atlas but excluded from the contrast, so it targets true ",
                  "cachexia; T/NK is dropped for too few cachectic cells).")
        ),
        hr(),
        h5("Current dataset"),
        tableOutput(ns("summary")),
        hr(),
        h5("Methods"),
        tags$ul(
          tags$li("Only raw (unfiltered) 10x matrices were released; an empty-droplet ",
                  "knee filter (min 200 genes / 500 UMIs, <20% mito) was applied per sample."),
          tags$li("Each dataset processed independently in scanpy: normalize + log1p, ",
                  "HVG → PCA → neighbours → Leiden → UMAP."),
          tags$li("Cell types annotated by marker-score argmax over canonical ",
                  "non-myofiber lineages; ", strong("myofiber/contaminant clusters were removed"),
                  " (this app shows only non-myofiber cells)."),
          tags$li("NF-κB reporter Neg/Pos fractions are pooled per condition."),
          tags$li(strong("Differential expression:"), " Tumor-vs-Control Wilcoxon within ",
                  "each cell type. ", strong("Marker genes:"), " one-vs-rest Wilcoxon."),
          tags$li(strong("Pathway enrichment:"), " GSEA (fgsea) ranks each cell type's ",
                  "Tumor-vs-Control genes by the Wilcoxon z-statistic against MSigDB ",
                  "Hallmark / Reactome / GO:BP (BH-adjusted within collection × cell type); ",
                  strong("PROGENy"), " scores 14 pathway activities per cell, reported as the ",
                  "cachexia − control change (effect size; pseudobulk p-values are ",
                  "exploratory as the NF-κB reporter fractions are not independent replicates). ",
                  "In both, a positive value = up in cachexia (red)."),
          tags$li("Expression values shown are log-normalized.")
        ),
        hr(),
        p(class = "text-muted",
          tags$em("Cell-type labels are marker-based defaults for exploration and ",
                  "should be validated experimentally. "),
          "GEO: ",
          tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248800",
                 target = "_blank", "GSE248800"), ".")
      )
    )
  )
}

aboutServer <- function(id, active) {
  moduleServer(id, function(input, output, session) {
    output$summary <- renderTable({
      d <- DS[[active()]]
      data.frame(
        Property = c("Dataset", "Cells (non-myofiber)", "Cell types",
                     "Conditions", "Genes"),
        Value = c(d$label, format(nrow(d$meta), big.mark = ","),
                  length(d$cell_types), paste(d$conditions, collapse = ", "),
                  format(length(d$gene_list), big.mark = ","))
      )
    }, striped = TRUE, spacing = "xs", width = "100%")
  })
}
