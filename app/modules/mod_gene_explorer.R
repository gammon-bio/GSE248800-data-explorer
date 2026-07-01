# mod_gene_explorer.R — per-gene UMAP + violin/dot across cell types.
# Dataset-aware: reads DS[[active()]] and pulls one gene-row from that dataset's
# HDF5 store. UMAP and expression plot live in separate sub-tabs (no stacking).

geneExplorerUI <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      title = "Gene & options", width = 320,
      div(class = "gene-search",
          selectizeInput(ns("gene"), "Gene", choices = NULL, multiple = FALSE,
                         options = list(placeholder = "Type a gene symbol...",
                                        maxOptions = 50))),
      radioButtons(ns("plot_type"), "Expression plot:",
                   choices = c("Violin" = "violin", "Dot plot" = "dot"),
                   selected = "violin"),
      checkboxInput(ns("split"), "Split by condition (Tumor vs Control)", value = TRUE),
      hr(),
      tags$strong("Downloads"),
      div(class = "d-flex gap-2 flex-wrap",
          downloadButton(ns("dl_umap"), "UMAP (PDF)", class = "btn-sm download-btn"),
          downloadButton(ns("dl_expr"), "Expression (PDF)", class = "btn-sm download-btn"))
    ),
    navset_card_tab(
      nav_panel("UMAP",
        card_body(class = "plot-container",
                  uiOutput(ns("umap_caption")),
                  plotlyOutput(ns("umap_plot"), height = "560px"))),
      nav_panel("Expression by cell type",
        card_body(class = "plot-container",
                  plotOutput(ns("expr_plot"), height = "480px")))
    )
  )
}

geneExplorerServer <- function(id, active) {
  moduleServer(id, function(input, output, session) {

    ds <- reactive(DS[[active()]])

    # Repopulate the gene dropdown when the dataset changes; keep the current
    # gene if it also exists in the newly selected dataset.
    observeEvent(active(), {
      gl <- ds()$gene_list
      sel <- isolate(input$gene)
      keep <- if (!is.null(sel) && sel %in% gl) sel else character(0)
      updateSelectizeInput(session, "gene", choices = gl, server = TRUE, selected = keep)
    }, ignoreInit = FALSE)

    gene_df <- reactive({
      req(input$gene)
      vec <- read_gene(active(), input$gene)
      if (is.null(vec)) {
        showNotification(sprintf("'%s' not in %s.", input$gene, ds()$label), type = "warning")
        return(NULL)
      }
      m <- ds()$meta
      data.frame(UMAP1 = m$UMAP1, UMAP2 = m$UMAP2, expression = vec,
                 cell_type = m$cell_type, condition = m$condition,
                 stringsAsFactors = FALSE)
    })

    output$umap_caption <- renderUI({
      if (is.null(input$gene) || input$gene == "")
        div(class = "text-muted small", "Select a gene to colour the UMAP by expression.")
      else div(class = "small", strong(input$gene), sprintf(" — log-norm expression in %s", ds()$label))
    })

    umap_gg <- reactive({
      df <- gene_df(); req(df)
      umap_expression_plot(df, input$gene, split_by = if (isTRUE(input$split)) "condition" else NULL)
    })
    expr_gg <- reactive({
      df <- gene_df(); req(df)
      if (input$plot_type == "violin") gene_violin(df, input$gene, split = isTRUE(input$split))
      else gene_dotplot(df, input$gene)
    })

    output$umap_plot <- renderPlotly({
      req(input$gene, gene_df())
      p <- plotly::toWebGL(ggplotly(umap_gg(), tooltip = "text"))
      p$x$data <- lapply(p$x$data, function(tr) { tr$hoveron <- NULL; tr })
      p
    })
    output$expr_plot <- renderPlot({ req(gene_df()); expr_gg() })

    output$dl_umap <- ggsave_handler(
      basename = reactive(paste0("umap_", input$gene %||% "gene", "_", active())),
      plot_fun = function() umap_gg(), width = 9, height = 7)
    output$dl_expr <- ggsave_handler(
      basename = reactive(paste0("expr_", input$gene %||% "gene", "_", active())),
      plot_fun = function() expr_gg(), width = 9, height = 6)
  })
}
