# app.R — entry point. global.R is auto-sourced by Shiny; ui.R and server.R are
# sourced here so the app can be launched with shiny::runApp("app").
source("ui.R", local = TRUE)
source("server.R", local = TRUE)
shinyApp(ui = ui, server = server)
