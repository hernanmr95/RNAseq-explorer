# R/qc_filtering.R
# Filtrado de genes de baja expresión y PCA exploratorio previo al modelado.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

#' Genera un PCA de las muestras a partir de una transformación de varianza estabilizada
#'
#' @param dds objeto DESeqDataSet (sin necesidad de haber corrido DESeq() aún)
#' @param condition_col columna de metadata usada para colorear los puntos
#' @param out_path ruta donde guardar el PNG del PCA
#' @return invisible(NULL); efecto secundario: guarda el plot en disco
run_exploratory_pca <- function(dds, condition_col, out_path) {
  # DESeq2::vst() usa por defecto una submuestra (nsub, mínimo 1000 genes)
  # para ajustar la tendencia de dispersión de forma rápida en datasets
  # grandes. Con pocos genes (como el dataset sintético de este repo, o
  # paneles dirigidos pequeños) eso falla, y la propia función recomienda
  # usar varianceStabilizingTransformation() directamente: hace el ajuste
  # completo sin submuestreo, más lento en datasets enormes pero correcto
  # para cualquier tamaño.
  vst_data <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  pca_data <- DESeq2::plotPCA(vst_data, intgroup = condition_col, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))

  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = .data[[condition_col]])) +
    geom_point(size = 3.5, alpha = 0.9) +
    labs(
      title = "PCA de muestras (transformación de varianza estabilizada)",
      x = sprintf("PC1 (%d%% varianza)", percent_var[1]),
      y = sprintf("PC2 (%d%% varianza)", percent_var[2])
    ) +
    theme_minimal(base_size = 12)

  ggsave(out_path, plot = p, width = 6, height = 5, dpi = 150)
  invisible(NULL)
}
