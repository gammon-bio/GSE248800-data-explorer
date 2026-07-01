# mod_cell_composition.R — dataset overview: UMAP by cell type / by condition,
# cell-type proportions by condition, and a counts table. Each in its own sub-tab.

cellCompositionUI <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    nav_panel("UMAP — cell type",
      card_body(class = "plot-container", plotlyOutput(ns("umap_ct"), height = "600px"))),
    nav_panel("UMAP — condition",
      card_body(class = "plot-container", plotlyOutput(ns("umap_cond"), height = "600px"))),
    nav_panel("Proportions",
      card_body(class = "plot-container", plotOutput(ns("bar"), height = "480px")),
      card_footer(downloadButton(ns("dl_bar"), "Download (PDF)", class = "btn-sm"))),
    nav_panel("Counts",
      card_body(DTOutput(ns("counts"))))
  )
}

cellCompositionServer <- function(id, active) {
  moduleServer(id, function(input, output, session) {

    ds <- reactive(DS[[active()]])

    # Subsample for interactive UMAPs (deterministic per dataset for stability).
    umap_sub <- reactive({
      m <- ds()$meta
      if (nrow(m) > 25000) {
        set.seed(1); m[sort(sample(nrow(m), 25000)), ]
      } else m
    })

    output$umap_ct <- renderPlotly({
      p <- umap_categorical_plot(umap_sub(), "cell_type", ds()$palette, "Cell type")
      plotly::toWebGL(ggplotly(p, tooltip = "text"))
    })

    output$umap_cond <- renderPlotly({
      cols <- condition_colors(ds()$conditions)
      p <- umap_categorical_plot(umap_sub(), "condition", cols, "Condition")
      plotly::toWebGL(ggplotly(p, tooltip = "text"))
    })

    bar_gg <- reactive({
      counts <- ds()$meta %>% count(condition, cell_type, name = "n")
      comp_bar(counts, ds()$palette)
    })
    output$bar <- renderPlot(bar_gg())

    output$counts <- renderDT({
      wide <- ds()$summary %>%
        select(cell_type, condition, n_cells) %>%
        pivot_wider(names_from = condition, values_from = n_cells, values_fill = 0)
      cond_cols <- setdiff(names(wide), "cell_type")
      wide$Total <- rowSums(wide[cond_cols])
      wide <- wide %>% arrange(desc(Total)) %>% rename(`Cell type` = cell_type)
      datatable(wide, options = list(pageLength = 20, dom = "t"), rownames = FALSE)
    })

    output$dl_bar <- ggsave_handler(
      basename = reactive(paste0("composition_", active())),
      plot_fun = function() bar_gg())
  })
}
