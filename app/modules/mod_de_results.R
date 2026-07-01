# mod_de_results.R — Tumor-vs-Control differential expression per cell type.
# de_results is a named list (cell_type -> data frame with gene, avg_log2FC,
# p_val, p_val_adj, pct.1 [tumor], pct.2 [control]). Volcano + table in sub-tabs.

deResultsUI <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      title = "DE filters", width = 320,
      selectInput(ns("cell_type"), "Cell type:", choices = NULL),
      sliderInput(ns("logfc"), "|log2FC| threshold:", min = 0, max = 3,
                  value = 0.25, step = 0.05),
      sliderInput(ns("pval"), "-log10(adj p) threshold:", min = 0, max = 50,
                  value = 1.3, step = 0.1),
      hr(),
      div(class = "d-flex gap-2 flex-wrap",
          downloadButton(ns("dl_table"), "Table (CSV)", class = "btn-sm download-btn"),
          downloadButton(ns("dl_volcano"), "Volcano (PDF)", class = "btn-sm download-btn"))
    ),
    navset_card_tab(
      nav_panel("Volcano",
        card_body(class = "plot-container", plotlyOutput(ns("volcano"), height = "560px"))),
      nav_panel("DE table",
        card_body(DTOutput(ns("table"))))
    )
  )
}

deResultsServer <- function(id, active) {
  moduleServer(id, function(input, output, session) {

    ds <- reactive(DS[[active()]])

    observeEvent(active(), {
      cts <- names(ds()$de_results)
      sel <- isolate(input$cell_type)
      keep <- if (!is.null(sel) && sel %in% cts) sel else cts[1]
      updateSelectInput(session, "cell_type", choices = cts, selected = keep)
    }, ignoreInit = FALSE)

    labels <- reactive({
      tum <- ds()$tumor_label
      list(tum = tum, up = paste("Up in", tum), down = paste("Down in", tum))
    })

    de_data <- reactive({
      req(input$cell_type, input$cell_type %in% names(ds()$de_results))
      df <- ds()$de_results[[input$cell_type]]
      df$neg_log10_padj <- -log10(pmax(df$p_val_adj, 1e-300))
      lab <- labels()
      sig <- abs(df$avg_log2FC) > input$logfc & df$neg_log10_padj > input$pval
      df$significance <- "NS"
      df$significance[sig & df$avg_log2FC > 0] <- lab$up
      df$significance[sig & df$avg_log2FC < 0] <- lab$down
      df
    })

    volcano_gg <- reactive({
      lab <- labels()
      de_volcano(de_data(),
                 title = sprintf("%s — %s vs Control (%s)", input$cell_type, lab$tum, ds()$label),
                 lfc_thresh = input$logfc, p_thresh = input$pval,
                 up_label = lab$up, down_label = lab$down)
    })

    output$volcano <- renderPlotly({
      ggplotly(volcano_gg(), tooltip = c("text", "x", "y")) %>%
        layout(legend = list(orientation = "h", y = 1.1))
    })

    output$table <- renderDT({
      lab <- labels()
      df <- de_data() %>%
        filter(significance != "NS") %>%
        select(gene, avg_log2FC, pct.1, pct.2, p_val_adj, significance) %>%
        arrange(p_val_adj) %>%
        mutate(avg_log2FC = round(avg_log2FC, 3),
               pct.1 = round(pct.1, 1), pct.2 = round(pct.2, 1),
               p_val_adj = signif(p_val_adj, 3))
      datatable(df,
                colnames = c("Gene", "log2FC", paste0("% ", lab$tum), "% Control",
                             "Adj. p", "Direction"),
                options = list(pageLength = 15, scrollX = TRUE),
                selection = "single", rownames = FALSE)
    })

    output$dl_table <- csv_handler(
      basename = reactive(paste0("DE_", gsub("[^0-9A-Za-z]+", "_", input$cell_type), "_", active())),
      df_fun = function() de_data() %>% filter(significance != "NS") %>%
        select(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj, significance) %>%
        arrange(p_val_adj))
    output$dl_volcano <- ggsave_handler(
      basename = reactive(paste0("volcano_", gsub("[^0-9A-Za-z]+", "_", input$cell_type), "_", active())),
      plot_fun = function() volcano_gg())
  })
}
