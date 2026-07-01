# ui.R — bslib page_navbar shell.
# A dataset toggle (Pancreatic / Colorectal / PDAC) lives in the navbar and
# drives every tab. fillable = FALSE makes tab bodies flow and SCROLL rather
# than being squeezed to the viewport; each data tab uses navset_card_tab
# sub-tabs so plots are shown one at a time instead of stacked.

ui <- bslib::page_navbar(
  title = "GSE248800 — Muscle Non-myofiber scRNA-seq Explorer",
  id = "main_nav",
  fillable = FALSE,
  theme = bslib::bs_theme(
    version = 5, bootswatch = "flatly",
    primary = "#2C3E50", "navbar-bg" = "#2C3E50"
  ),
  header = tags$head(
    waiter::useWaiter(),
    includeCSS("www/custom.css")
  ),

  nav_panel("Gene Explorer",        icon = icon("dna"),          geneExplorerUI("gene")),
  nav_panel("Differential Expression", icon = icon("chart-bar"), deResultsUI("de")),
  nav_panel("Pathway Enrichment",   icon = icon("diagram-project"), pathwaysUI("path")),
  nav_panel("Cell Composition",     icon = icon("circle-nodes"), cellCompositionUI("comp")),
  nav_panel("Marker Genes",         icon = icon("list"),         markerGenesUI("markers")),
  nav_panel("About",                icon = icon("info-circle"),  aboutUI("about")),

  # --- Dataset toggle pinned to the right of the navbar ---------------------
  nav_spacer(),
  nav_item(
    div(
      class = "dataset-toggle",
      shinyWidgets::radioGroupButtons(
        inputId = "active_dataset",
        label = NULL,
        choiceNames = unname(DATASET_LABELS[AVAILABLE]),
        choiceValues = AVAILABLE,
        selected = DEFAULT_DATASET,
        size = "sm",
        status = "light"
      )
    )
  )
)
