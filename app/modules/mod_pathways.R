# mod_pathways.R — Pathway Enrichment tab.
# Displays precomputed GSEA (fgsea: Hallmark / Reactome / GO:BP) and PROGENy
# pathway-activity results per cell type for the Tumor-vs-Control (cachexia)
# contrast. Everything is precomputed offline in data-raw/04_pathways.R; this
# module only filters tidy data frames and draws diverging heatmaps/bars.
# Sign convention everywhere: value > 0 = UP in cachexia = red.

COLL_CHOICES <- c("Hallmark (H)" = "H",
                  "Reactome (CP)" = "C2_CP_REACTOME",
                  "GO Biological Process" = "C5_GO_BP")

pathwaysUI <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      title = "Pathway filters", width = 320,
      selectInput(ns("collection"), "Gene-set collection:", choices = COLL_CHOICES),
      selectInput(ns("cell_type"), "Cell type:", choices = NULL),
      sliderInput(ns("top_n"), "Top N pathways / direction:", min = 5, max = 30, value = 10),
      sliderInput(ns("padj_cut"), "padj cutoff (significance *):", min = 0, max = 0.25,
                  value = 0.05, step = 0.01),
      checkboxInput(ns("sig_only"), "Significant only", value = FALSE),
      hr(),
      selectizeInput(ns("focus_pathway"), "Pathway (cross-cell-type view):",
                     choices = NULL, options = list(placeholder = "Type a pathway...")),
      radioButtons(ns("result_type"), "Table shows:",
                   choices = c("GSEA", "PROGENy"), inline = TRUE),
      hr(),
      div(class = "d-flex gap-2 flex-wrap",
          downloadButton(ns("dl_heatmap"), "Heatmap (PDF)", class = "btn-sm download-btn"),
          downloadButton(ns("dl_bar"), "Bars (PDF)", class = "btn-sm download-btn"),
          downloadButton(ns("dl_progeny"), "PROGENy (PDF)", class = "btn-sm download-btn"),
          downloadButton(ns("dl_csv"), "Table (CSV)", class = "btn-sm download-btn"))
    ),
    navset_card_tab(
      nav_panel("NES heatmap",
        card_body(class = "plot-container",
                  div(class = "small text-muted", "Top enriched pathways per cell type (red = up in cachexia). * = padj < cutoff."),
                  plotlyOutput(ns("heatmap"), height = "600px"))),
      nav_panel("Top pathways",
        card_body(class = "plot-container",
                  plotlyOutput(ns("bar"), height = "600px"))),
      nav_panel("Pathway across cell types",
        card_body(class = "plot-container",
                  plotlyOutput(ns("crosscut"), height = "520px"))),
      nav_panel("PROGENy activity",
        card_body(class = "plot-container",
                  div(class = "small text-muted", "PROGENy 14-pathway activity change (cachexia − control). * = pseudobulk padj < 0.05 (exploratory)."),
                  plotlyOutput(ns("progeny"), height = "560px"))),
      nav_panel("Table",
        card_body(DTOutput(ns("table"))))
    )
  )
}

pathwaysServer <- function(id, active) {
  moduleServer(id, function(input, output, session) {

    ds     <- reactive(DS[[active()]])
    gsea   <- reactive(ds()$gsea)       # may be NULL on a partial build
    prog   <- reactive(ds()$progeny)
    tumor  <- reactive(ds()$tumor_label)

    # Repopulate cell-type + pathway choices when the dataset/collection changes.
    observeEvent(list(active(), input$collection), {
      g <- gsea()
      cts <- if (is.null(g)) ds()$cell_types else sort(unique(as.character(g$cell_type)))
      sel_ct <- isolate(input$cell_type)
      updateSelectInput(session, "cell_type", choices = cts,
                        selected = if (!is.null(sel_ct) && sel_ct %in% cts) sel_ct else cts[1])
      paths <- if (is.null(g)) character(0)
               else sort(unique(g$pathway_label[as.character(g$collection) == input$collection]))
      sel_p <- isolate(input$focus_pathway)
      keep_p <- if (!is.null(sel_p) && sel_p %in% paths) sel_p else paths[1]
      updateSelectizeInput(session, "focus_pathway", choices = paths, selected = keep_p, server = TRUE)
    }, ignoreInit = FALSE)

    # GSEA rows for the active collection (character-cast for safe filtering).
    gsea_coll <- reactive({
      g <- gsea(); validate(need(!is.null(g) && nrow(g) > 0,
        "Pathway enrichment has not been precomputed for this dataset."))
      g <- g[as.character(g$collection) == input$collection, , drop = FALSE]
      g <- g[is.finite(g$NES), , drop = FALSE]   # fgsea emits NA NES for a few sets; drop them
      g$cell_type <- as.character(g$cell_type)
      g
    })

    gsea_tooltip <- function(d) paste0(
      "<b>", d$pathway_label, "</b><br>", d$cell_type,
      "<br>NES: ", sprintf("%.2f", d$NES), "<br>padj: ", sprintf("%.2e", d$padj),
      "<br>size: ", d$size, "<br>lead: ", substr(d$leadingEdge, 1, 60))

    # Top-N up + down per cell type -> the columns shown in the overview heatmap.
    heatmap_df <- reactive({
      g <- gsea_coll(); n <- input$top_n
      cand <- if (isTRUE(input$sig_only)) g[!is.na(g$padj) & g$padj < input$padj_cut, ] else g
      validate(need(nrow(cand) > 0, "No pathways pass the current filters."))
      picks <- unique(unlist(lapply(split(cand, cand$cell_type), function(d) {
        up <- utils::head(d$pathway[order(-d$NES)], n)
        dn <- utils::head(d$pathway[order(d$NES)], n)
        c(up, dn)
      })))
      d <- g[g$pathway %in% picks, , drop = FALSE]
      # Cap columns for readability: the big collections (Reactome, GO:BP) can
      # union to 100+ pathways. Keep the most discriminative ones (largest
      # |NES| across cell types) so the overview stays legible; the Top-pathways
      # and Table sub-tabs give the full per-cell-type list.
      MAX_COLS <- 40
      if (length(unique(d$pathway_label)) > MAX_COLS) {
        keep_pw <- names(utils::head(sort(tapply(abs(d$NES), d$pathway_label,
                                                 max, na.rm = TRUE), decreasing = TRUE), MAX_COLS))
        d <- d[d$pathway_label %in% keep_pw, , drop = FALSE]
      }
      # order columns + rows by mean NES so the inflammatory (red) block clusters
      pw_order <- names(sort(tapply(d$NES, d$pathway_label, mean, na.rm = TRUE)))
      ct_order <- names(sort(tapply(d$NES, d$cell_type, mean, na.rm = TRUE)))
      # Complete the grid (every shown pathway x every cell type) so ggplotly
      # gets a full rectangle; genuinely-absent cells render grey (NA).
      full <- expand.grid(pathway_label = pw_order, cell_type = ct_order,
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      d <- merge(full, d[, c("pathway_label", "cell_type", "NES", "padj", "size", "leadingEdge")],
                 by = c("pathway_label", "cell_type"), all.x = TRUE)
      data.frame(
        xvar = factor(d$pathway_label, levels = pw_order),
        yvar = factor(d$cell_type, levels = ct_order),
        value = d$NES,
        sig = !is.na(d$padj) & d$padj < input$padj_cut,
        tooltip = gsea_tooltip(d), stringsAsFactors = FALSE)
    })

    heatmap_gg <- reactive(enrich_heatmap(heatmap_df(), "NES"))
    output$heatmap <- renderPlotly({
      ggplotly(heatmap_gg(), tooltip = "text") %>% layout(margin = list(b = 120))
    })

    # Top up/down bars for the selected cell type.
    bar_df <- reactive({
      req(input$cell_type)
      g <- gsea_coll(); d <- g[g$cell_type == input$cell_type, , drop = FALSE]
      if (isTRUE(input$sig_only)) d <- d[!is.na(d$padj) & d$padj < input$padj_cut, ]
      validate(need(nrow(d) > 0, "No pathways for this cell type / filters."))
      n <- input$top_n
      d <- rbind(utils::head(d[order(-d$NES), ], n), utils::head(d[order(d$NES), ], n))
      d <- d[!duplicated(d$pathway), ]
      data.frame(label = factor(d$pathway_label, levels = d$pathway_label[order(d$NES)]),
                 value = d$NES, tooltip = gsea_tooltip(d), stringsAsFactors = FALSE)
    })
    bar_gg <- reactive(pathway_bar(bar_df(), "NES",
      sprintf("%s — %s vs Control", input$cell_type, tumor())))
    output$bar <- renderPlotly(ggplotly(bar_gg(), tooltip = "text"))

    # One pathway across all cell types.
    crosscut_df <- reactive({
      req(input$focus_pathway)
      g <- gsea_coll(); d <- g[g$pathway_label == input$focus_pathway, , drop = FALSE]
      validate(need(nrow(d) > 0, "Select a pathway."))
      data.frame(label = factor(d$cell_type, levels = d$cell_type[order(d$NES)]),
                 value = d$NES, tooltip = gsea_tooltip(d), stringsAsFactors = FALSE)
    })
    crosscut_gg <- reactive(pathway_bar(crosscut_df(), "NES", input$focus_pathway))
    output$crosscut <- renderPlotly(ggplotly(crosscut_gg(), tooltip = "text"))

    # PROGENy activity heatmap (cell type x 14 pathways).
    progeny_df <- reactive({
      p <- prog(); validate(need(!is.null(p) && nrow(p) > 0,
        "PROGENy has not been precomputed for this dataset."))
      data.frame(
        xvar = factor(as.character(p$pathway)),
        yvar = factor(as.character(p$cell_type),
                      levels = names(sort(tapply(p$contrast, as.character(p$cell_type), mean, na.rm = TRUE)))),
        value = p$contrast,
        sig = !is.na(p$padj) & p$padj < 0.05 & !p$low_n,
        tooltip = paste0("<b>", p$pathway, "</b><br>", p$cell_type,
                         "<br>contrast: ", sprintf("%+.2f", p$contrast),
                         "<br>Cohen's d: ", sprintf("%.2f", p$cohens_d),
                         "<br>padj: ", ifelse(is.na(p$padj), "NA", sprintf("%.2e", p$padj))),
        stringsAsFactors = FALSE)
    })
    progeny_gg <- reactive(enrich_heatmap(progeny_df(), "Activity\nΔ (cachexia)"))
    output$progeny <- renderPlotly({
      ggplotly(progeny_gg(), tooltip = "text") %>% layout(margin = list(b = 80))
    })

    # Table (GSEA or PROGENy).
    output$table <- renderDT({
      if (input$result_type == "GSEA") {
        g <- gsea_coll(); d <- g[g$cell_type == input$cell_type, , drop = FALSE]
        d <- d[order(d$padj), c("pathway_label", "NES", "padj", "size", "direction", "leadingEdge")]
        d$NES <- round(d$NES, 2); d$padj <- signif(d$padj, 3)
        datatable(d, colnames = c("Pathway", "NES", "padj", "Size", "Direction", "Leading edge"),
                  options = list(pageLength = 15, scrollX = TRUE), selection = "single", rownames = FALSE)
      } else {
        p <- prog(); validate(need(!is.null(p), "PROGENy not precomputed."))
        d <- p[order(-abs(p$contrast)), c("cell_type", "pathway", "contrast", "cohens_d", "padj", "low_n")]
        d$contrast <- round(d$contrast, 3); d$cohens_d <- round(d$cohens_d, 3); d$padj <- signif(d$padj, 3)
        datatable(d, colnames = c("Cell type", "Pathway", "Contrast", "Cohen's d", "padj (pseudobulk)", "low n"),
                  options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
      }
    })

    # Selecting a GSEA table row focuses that pathway in the cross-cut tab.
    observeEvent(input$table_rows_selected, {
      if (input$result_type != "GSEA") return()
      g <- gsea_coll(); d <- g[g$cell_type == input$cell_type, , drop = FALSE]
      d <- d[order(d$padj), ]
      pw <- d$pathway_label[input$table_rows_selected]
      if (length(pw)) updateSelectizeInput(session, "focus_pathway", selected = pw)
    })

    # ---- downloads ----
    base <- reactive(paste0(input$collection, "_", gsub("[^0-9A-Za-z]+", "_", input$cell_type %||% ""), "_", active()))
    output$dl_heatmap <- ggsave_handler(reactive(paste0("gsea_heatmap_", base())),
                                        function() heatmap_gg(), width = 13, height = 7)
    output$dl_bar     <- ggsave_handler(reactive(paste0("gsea_bars_", base())),
                                        function() bar_gg(), width = 10, height = 8)
    output$dl_progeny <- ggsave_handler(reactive(paste0("progeny_", active())),
                                        function() progeny_gg(), width = 9, height = 6)
    output$dl_csv <- csv_handler(reactive(paste0("pathways_", input$result_type, "_", active())),
      function() if (input$result_type == "GSEA") gsea_coll() else prog())
  })
}
