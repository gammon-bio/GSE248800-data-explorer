# plot_helpers.R — PURE ggplot builders shared by render + download.
# Inputs are plain data frames (already aggregated / joined by the module); no
# reactivity, no globals. Continuous expression uses magma viridis; Control is
# blue (#3498DB) and tumour conditions are red/orange, matching zhang_shiny.

CTRL_COLOR  <- "#3498DB"
TUMOR_COLOR <- "#E74C3C"

# condition_colors(levels) — Control -> blue, other conditions -> a red/orange
# ramp (handles the human 3-level Control / PDAC_weight_stable / PDAC_cachectic).
condition_colors <- function(levels) {
  levels <- as.character(levels)
  others <- setdiff(levels, "Control")
  ramp <- if (length(others) <= 1) TUMOR_COLOR
          else grDevices::colorRampPalette(c("#E74C3C", "#F39C12"))(length(others))
  stats::setNames(c(CTRL_COLOR, ramp), c("Control", others))[levels]
}

# UMAP coloured by continuous gene expression (high expressers drawn on top).
umap_expression_plot <- function(df, gene_label, split_by = NULL) {
  df <- df[order(df$expression), , drop = FALSE]
  df$tooltip <- paste0("Cell type: ", df$cell_type,
                       "<br>", gene_label, ": ", round(df$expression, 2),
                       "<br>Condition: ", df$condition)
  p <- ggplot2::ggplot(df, ggplot2::aes(UMAP1, UMAP2, color = expression, text = tooltip)) +
    ggplot2::geom_point(size = 0.4, alpha = 0.85) +
    viridis::scale_color_viridis(option = "magma", name = gene_label) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right")
  if (!is.null(split_by)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~", split_by)))
  p
}

# UMAP coloured by a categorical column (cell_type or condition).
umap_categorical_plot <- function(df, color_col, colors, legend_title) {
  df$.col <- as.character(df[[color_col]])
  df$tooltip <- paste0(legend_title, ": ", df$.col)
  ggplot2::ggplot(df, ggplot2::aes(UMAP1, UMAP2, color = .col, text = tooltip)) +
    ggplot2::geom_point(size = 0.35, alpha = 0.8) +
    ggplot2::scale_color_manual(values = colors, name = legend_title) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
}

# Violin of one gene across cell types, optionally split by condition.
gene_violin <- function(df, gene_label, split = TRUE) {
  cond_cols <- condition_colors(sort(unique(df$condition)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = cell_type, y = expression))
  if (split) {
    p <- p + ggplot2::geom_violin(ggplot2::aes(fill = condition),
                                  scale = "width", position = ggplot2::position_dodge(0.9),
                                  linewidth = 0.2) +
      ggplot2::scale_fill_manual(values = cond_cols, name = "Condition")
  } else {
    p <- p + ggplot2::geom_violin(fill = "#7f8c8d", scale = "width", linewidth = 0.2)
  }
  p +
    ggplot2::labs(x = NULL, y = sprintf("%s (log-norm)", gene_label)) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "top")
}

# Dot plot: mean expression + % expressing per cell type (one gene).
gene_dotplot <- function(df, gene_label) {
  agg <- dplyr::summarise(dplyr::group_by(df, cell_type),
                          avg_exp = mean(expression),
                          pct_exp = mean(expression > 0) * 100, .groups = "drop")
  agg$cell_type <- factor(agg$cell_type, levels = agg$cell_type[order(agg$avg_exp)])
  agg$gene <- gene_label
  ggplot2::ggplot(agg, ggplot2::aes(gene, cell_type)) +
    ggplot2::geom_point(ggplot2::aes(size = pct_exp, color = avg_exp)) +
    viridis::scale_color_viridis(option = "magma", name = "Mean expr") +
    ggplot2::scale_size(name = "% expressing", range = c(1, 11)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "right")
}

# Volcano for a Tumor-vs-Control DE table (columns: gene, avg_log2FC,
# neg_log10_padj, significance). Colours: up=red, down=blue, NS=grey.
de_volcano <- function(df, title, lfc_thresh, p_thresh, up_label, down_label) {
  sig_colors <- stats::setNames(c(TUMOR_COLOR, CTRL_COLOR, "#BDC3C7"),
                                c(up_label, down_label, "NS"))
  ggplot2::ggplot(df, ggplot2::aes(avg_log2FC, neg_log10_padj,
                                   color = significance, text = gene)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.4) +
    ggplot2::scale_color_manual(values = sig_colors, name = NULL) +
    ggplot2::geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
                        linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = p_thresh, linetype = "dashed", color = "grey50") +
    ggplot2::labs(x = "log2 fold change", y = "-log10(adjusted p-value)", title = title) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "top")
}

# Stacked proportion bar of cell-type composition by condition.
comp_bar <- function(counts_df, palette) {
  ggplot2::ggplot(counts_df, ggplot2::aes(x = condition, y = n, fill = cell_type)) +
    ggplot2::geom_col(position = "fill", width = 0.65) +
    ggplot2::scale_fill_manual(values = palette, name = "Cell type") +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(x = NULL, y = "Proportion") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "right")
}

# Heatmap of per-gene scaled mean expression across cell types.
marker_heatmap <- function(heat_df) {
  ggplot2::ggplot(heat_df, ggplot2::aes(gene, cell_type, fill = scaled_expr)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    viridis::scale_fill_viridis(option = "magma", name = "Scaled\nexpr") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
                   axis.text.y = ggplot2::element_text(size = 10),
                   panel.grid = ggplot2::element_blank())
}

# ---- Pathway Enrichment builders -------------------------------------------
# Signed quantities (NES, PROGENy contrast) need a DIVERGING scale: Control/down
# = blue, tumour/up = red, 0 = neutral grey (matches the volcano's "red = up in
# cachexia"). Magma stays reserved for non-negative expression.
diverging_fill <- function(name, limits = NULL) {
  ggplot2::scale_fill_gradient2(
    low = CTRL_COLOR, mid = "grey95", high = TUMOR_COLOR, midpoint = 0,
    name = name, limits = limits, oob = scales::squish, na.value = "grey85")
}

# Heatmap for enrichment. df columns: xvar (factor), yvar (factor), value (num),
# sig (lgl), tooltip (chr). Used for both the GSEA NES heatmap and PROGENy.
enrich_heatmap <- function(df, fill_name, xlab = NULL, ylab = NULL) {
  lim <- stats::quantile(abs(df$value), 0.98, na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- max(abs(df$value), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  p <- ggplot2::ggplot(df, ggplot2::aes(xvar, yvar, fill = value, text = tooltip)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    diverging_fill(fill_name, limits = c(-lim, lim)) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                   axis.text.y = ggplot2::element_text(size = 10),
                   panel.grid = ggplot2::element_blank())
  sigdf <- df[!is.na(df$sig) & df$sig, , drop = FALSE]
  if (nrow(sigdf))
    p <- p + ggplot2::geom_text(data = sigdf, ggplot2::aes(xvar, yvar, label = "*"),
                                inherit.aes = FALSE, size = 4.5, vjust = 0.72)
  p
}

# Horizontal bar of a signed statistic. df columns: label (factor, pre-ordered),
# value (num), tooltip (chr). Used for top up/down pathways and pathway-across-
# cell-types.
pathway_bar <- function(df, xlab = "NES", title = NULL) {
  lim <- max(abs(df$value), na.rm = TRUE); if (!is.finite(lim) || lim == 0) lim <- 1
  ggplot2::ggplot(df, ggplot2::aes(value, label, fill = value, text = tooltip)) +
    ggplot2::geom_col() +
    diverging_fill(xlab, limits = c(-lim, lim)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey40") +
    ggplot2::labs(x = xlab, y = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 9))
}

# Dot plot of a custom gene set across cell types.
marker_dotplot <- function(dot_df) {
  ggplot2::ggplot(dot_df, ggplot2::aes(gene, cell_type)) +
    ggplot2::geom_point(ggplot2::aes(size = pct_exp, color = avg_exp)) +
    viridis::scale_color_viridis(option = "magma", name = "Avg\nexpr") +
    ggplot2::scale_size_continuous(name = "% expressing", range = c(1, 10)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 11),
                   axis.text.y = ggplot2::element_text(size = 10),
                   panel.grid.major = ggplot2::element_line(color = "grey90"))
}
