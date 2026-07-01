# mod_marker_genes.R — top one-vs-rest marker heatmap + custom-gene dot plot.
# Expression is aggregated per cell type on demand from the dataset's HDF5 store
# (read_genes). Heatmap and dot plot live in separate sub-tabs.

markerGenesUI <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      title = "Options", width = 320,
      sliderInput(ns("n_markers"), "Top N markers per cell type:",
                  min = 3, max = 15, value = 5, step = 1),
      hr(),
      h6("Custom gene set (dot plot)"),
      selectizeInput(ns("custom_genes"), "Genes:", choices = NULL, multiple = TRUE,
                     options = list(placeholder = "Type gene symbols...",
                                    maxOptions = 50, maxItems = 30)),
      actionButton(ns("load_defaults"), "Load lineage markers",
                   class = "btn-sm btn-outline-primary"),
      hr(),
      div(class = "d-flex gap-2 flex-wrap",
          downloadButton(ns("dl_heatmap"), "Heatmap (PDF)", class = "btn-sm download-btn"),
          downloadButton(ns("dl_dot"), "Dot plot (PDF)", class = "btn-sm download-btn"))
    ),
    navset_card_tab(
      nav_panel("Marker heatmap",
        card_body(class = "plot-container", plotOutput(ns("heatmap"), height = "560px"))),
      nav_panel("Custom dot plot",
        card_body(class = "plot-container", plotOutput(ns("dotplot"), height = "560px")))
    )
  )
}

# Canonical lineage markers used by "Load lineage markers" (case-insensitively
# matched to the active dataset's gene list, so it works for mouse and human).
.DEFAULT_MARKERS <- c("Ptprc", "Adgre1", "Cd68", "Csf1r", "Cd3e", "Nkg7", "Ms4a1",
                      "S100a8", "Pdgfra", "Dcn", "Col1a1", "Pecam1", "Cdh5",
                      "Rgs5", "Acta2", "Pax7", "Adipoq", "Plp1")

markerGenesServer <- function(id, active) {
  moduleServer(id, function(input, output, session) {

    ds <- reactive(DS[[active()]])
    cell_type_levels <- reactive(ds()$cell_types)

    observeEvent(active(), {
      updateSelectizeInput(session, "custom_genes", choices = ds()$gene_list,
                           server = TRUE, selected = isolate(input$custom_genes))
    }, ignoreInit = FALSE)

    observeEvent(input$load_defaults, {
      gl <- ds()$gene_list
      hits <- gl[toupper(gl) %in% toupper(.DEFAULT_MARKERS)]
      updateSelectizeInput(session, "custom_genes", selected = hits)
    })

    # Per-cell-type mean of a set of genes -> matrix (cell_types x genes).
    ct_means <- function(genes) {
      exprs <- read_genes(active(), genes)
      genes <- names(exprs)
      if (!length(genes)) return(NULL)
      ct <- factor(ds()$meta$cell_type, levels = cell_type_levels())
      mat <- vapply(genes, function(g) tapply(exprs[[g]], ct, mean), numeric(nlevels(ct)))
      list(mat = mat, genes = genes, exprs = exprs, ct = ct)
    }

    top_markers <- reactive({
      ds()$marker_genes %>% group_by(cluster) %>%
        slice_max(order_by = avg_log2FC, n = input$n_markers, with_ties = FALSE) %>%
        ungroup()
    })

    heatmap_gg <- reactive({
      td <- top_markers()
      genes <- unique(td$gene)
      cm <- ct_means(genes); req(cm)
      scaled <- scale(cm$mat)                 # z-score per gene (column)
      scaled[is.nan(scaled)] <- 0
      gene_order <- td %>% arrange(cluster) %>% pull(gene) %>% unique()
      gene_order <- gene_order[gene_order %in% cm$genes]
      heat_df <- as.data.frame(scaled) %>%
        mutate(cell_type = rownames(scaled)) %>%
        pivot_longer(-cell_type, names_to = "gene", values_to = "scaled_expr")
      heat_df$gene <- factor(heat_df$gene, levels = gene_order)
      heat_df$cell_type <- factor(heat_df$cell_type, levels = cell_type_levels())
      marker_heatmap(heat_df)
    })
    output$heatmap <- renderPlot({ p <- heatmap_gg(); if (!is.null(p)) p })

    dotplot_gg <- reactive({
      req(input$custom_genes)
      cm <- ct_means(input$custom_genes); req(cm)
      pct <- vapply(cm$genes, function(g) tapply(cm$exprs[[g]] > 0, cm$ct, mean) * 100,
                    numeric(nlevels(cm$ct)))
      dot_df <- data.frame(
        cell_type = rep(rownames(cm$mat), times = length(cm$genes)),
        gene = rep(cm$genes, each = nrow(cm$mat)),
        avg_exp = as.vector(cm$mat), pct_exp = as.vector(pct))
      dot_df$gene <- factor(dot_df$gene, levels = cm$genes)
      dot_df$cell_type <- factor(dot_df$cell_type, levels = rev(cell_type_levels()))
      marker_dotplot(dot_df)
    })
    output$dotplot <- renderPlot(dotplot_gg())

    output$dl_heatmap <- ggsave_handler(
      basename = reactive(paste0("marker_heatmap_", active())),
      plot_fun = function() heatmap_gg(), width = 12, height = 7)
    output$dl_dot <- ggsave_handler(
      basename = reactive(paste0("marker_dotplot_", active())),
      plot_fun = function() dotplot_gg(), width = 10, height = 7)
  })
}
