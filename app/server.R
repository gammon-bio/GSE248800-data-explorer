# server.R — resolves the active dataset from the navbar toggle and hands a
# reactive dataset key to every feature module. Modules pull their data from the
# global DS registry (DS[[key]]) and read expression via read_gene(key, gene).

server <- function(input, output, session) {

  # Active dataset key (falls back to the default if unset / invalid).
  active <- reactive({
    k <- input$active_dataset
    if (is.null(k) || !(k %in% names(DS))) DEFAULT_DATASET else k
  })

  geneExplorerServer("gene", active)
  deResultsServer("de", active)
  pathwaysServer("path", active)
  cellCompositionServer("comp", active)
  markerGenesServer("markers", active)
  aboutServer("about", active)
}
