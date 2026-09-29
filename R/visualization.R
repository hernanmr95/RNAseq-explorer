# R/visualization.R
# Genera las visualizaciones estándar de un análisis de expresión diferencial:
# volcano plot y heatmap de los genes más variables.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

#' Volcano plot con los genes significativos resaltados y los top N etiquetados
#'
#' @param results_df tibble/data.frame con columnas log2FoldChange, padj, status, gene_id
#' @param out_path ruta del PNG de salida
#' @param top_n número de genes a etiquetar por nombre (los de menor padj)
#' @return invisible(NULL)
plot_volcano <- function(results_df, out_path, top_n = 8) {
  results_df <- results_df[!is.na(results_df$padj), ]
  top_genes <- head(results_df[order(results_df$padj), ], top_n)

  p <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
    geom_point(alpha = 0.7, size = 1.8) +
    scale_color_manual(values = c(UP = "#d73027", DOWN = "#4575b4", NS = "grey70")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(
      data = top_genes, aes(label = gene_id),
      size = 3, max.overlaps = 20, show.legend = FALSE
    ) +
    labs(
      title = "Expresión diferencial: tratado vs. control",
      x = "log2 Fold Change", y = "-log10(padj)", color = "Estado"
    ) +
    theme_minimal(base_size = 12)

  ggsave(out_path, plot = p, width = 7, height = 6, dpi = 150)
  invisible(NULL)
}

#' Heatmap de los genes más variables (sobre valores normalizados/VST)
#'
#' @param vst_matrix matriz de expresión transformada (genes x muestras)
#' @param metadata data.frame de metadata (para la anotación de columnas)
#' @param condition_col columna usada como anotación de color
#' @param out_path ruta del PNG de salida
#' @param n_genes número de genes más variables a mostrar
#' @return invisible(NULL)
plot_top_variable_heatmap <- function(vst_matrix, metadata, condition_col, out_path, n_genes = 20) {
  gene_variance <- apply(vst_matrix, 1, var)
  top_genes <- names(sort(gene_variance, decreasing = TRUE))[seq_len(min(n_genes, length(gene_variance)))]

  mat <- vst_matrix[top_genes, , drop = FALSE]
  mat_scaled <- t(scale(t(mat)))  # z-score por gen, para comparar patrones entre muestras

  annotation_col <- data.frame(condition = metadata[[condition_col]])
  rownames(annotation_col) <- metadata$sample_id

  pheatmap::pheatmap(
    mat_scaled,
    annotation_col = annotation_col,
    show_rownames = TRUE,
    fontsize_row = 7,
    filename = out_path,
    width = 6.5,
    height = 6
  )
  invisible(NULL)
}
