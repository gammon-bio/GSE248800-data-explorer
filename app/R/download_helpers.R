# download_helpers.R — shared download handlers (from kpc_tme).
# plot_fun / df_fun are zero-arg functions so renderPlot/renderDT and the
# download share one source of truth.

stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

# basename may be a plain string OR a function/reactive (evaluated at download
# time) so filenames can reflect the current gene / dataset.
.resolve_base <- function(basename) if (is.function(basename)) basename() else basename

ggsave_handler <- function(basename, plot_fun, device = "pdf",
                           width = 9, height = 6, dpi = 150) {
  shiny::downloadHandler(
    filename = function() sprintf("%s_%s.%s", .resolve_base(basename), stamp(), device),
    content  = function(file) {
      p <- plot_fun()
      if (!is.null(p)) {
        ggplot2::ggsave(file, plot = p, device = device,
                        width = width, height = height, dpi = dpi, units = "in")
      }
    }
  )
}

csv_handler <- function(basename, df_fun) {
  shiny::downloadHandler(
    filename = function() sprintf("%s_%s.csv", .resolve_base(basename), stamp()),
    content  = function(file) utils::write.csv(df_fun(), file, row.names = FALSE)
  )
}

# Null-coalescing helper (guard if a module is sourced standalone).
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a
}
