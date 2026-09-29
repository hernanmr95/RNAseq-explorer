# R/deseq2_analysis.R
# Ajuste del modelo DESeq2 y extracción de la tabla de resultados del contraste
# de interés (condición tratada vs. referencia).

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

#' Construye el DESeqDataSet a partir de la matriz de counts filtrada y la metadata
#'
#' @param count_matrix matriz de conteos (genes x muestras), ya filtrada
#' @param metadata data.frame con las columnas de diseño experimental
#' @param condition_col nombre de la columna que define los grupos a comparar
#' @param reference_level nivel de referencia (p. ej. "control")
#' @return objeto DESeqDataSet listo para DESeq2::DESeq()
build_dds <- function(count_matrix, metadata, condition_col, reference_level) {
  metadata[[condition_col]] <- relevel(
    factor(metadata[[condition_col]]),
    ref = reference_level
  )

  design_formula <- stats::as.formula(paste("~", condition_col))

  DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData   = metadata,
    design    = design_formula
  )
}

#' Ejecuta DESeq2 y devuelve la tabla de resultados como tibble ordenado por padj
#'
#' @param dds DESeqDataSet (salida de build_dds)
#' @param condition_col columna de condición usada en el diseño
#' @param reference_level nivel de referencia usado en el contraste
#' @return tibble con gene_id, baseMean, log2FoldChange, pvalue, padj
run_deseq2 <- function(dds, condition_col, reference_level) {
  dds <- DESeq(dds, quiet = TRUE)

  treatment_level <- setdiff(levels(dds[[condition_col]]), reference_level)[1]
  contrast <- c(condition_col, treatment_level, reference_level)

  res <- results(dds, contrast = contrast)

  as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    arrange(padj)
}
