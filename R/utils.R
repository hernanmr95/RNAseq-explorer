# R/utils.R
# Funciones auxiliares compartidas por los distintos pasos del pipeline.

#' Valida que las muestras de la matriz de counts coincidan con las de la metadata
#'
#' @param counts data.frame de conteos, con genes en filas (primera columna = gene_id)
#' @param metadata data.frame con al menos las columnas sample_id y la de condición
#' @return TRUE invisible si todo es correcto; lanza error en caso contrario
validate_inputs <- function(counts, metadata) {
  sample_cols <- colnames(counts)[-1]
  meta_samples <- metadata$sample_id

  missing_in_meta <- setdiff(sample_cols, meta_samples)
  missing_in_counts <- setdiff(meta_samples, sample_cols)

  if (length(missing_in_meta) > 0) {
    stop(sprintf(
      "Las siguientes muestras están en la matriz de counts pero no en la metadata: %s",
      paste(missing_in_meta, collapse = ", ")
    ))
  }
  if (length(missing_in_counts) > 0) {
    stop(sprintf(
      "Las siguientes muestras están en la metadata pero no en la matriz de counts: %s",
      paste(missing_in_counts, collapse = ", ")
    ))
  }
  invisible(TRUE)
}

#' Filtra genes de baja expresión
#'
#' Mantiene únicamente los genes cuya suma de counts across todas las
#' muestras supera un umbral mínimo. Es un filtrado conservador previo
#' al análisis con DESeq2 (que ya hace su propio filtrado interno, pero
#' reducir la matriz de entrada acelera el ajuste del modelo).
#'
#' @param count_matrix matriz numérica de conteos (genes x muestras)
#' @param min_total_count numérico, suma mínima de counts por gen
#' @return matriz filtrada
filter_low_expression <- function(count_matrix, min_total_count = 10) {
  keep <- rowSums(count_matrix) >= min_total_count
  count_matrix[keep, , drop = FALSE]
}

#' Etiqueta genes como up/down/no significativos según log2FC y padj
#'
#' @param results_df data.frame con columnas log2FoldChange y padj
#' @param lfc_threshold umbral absoluto de log2 fold-change
#' @param padj_threshold umbral de p-valor ajustado
#' @return el mismo data.frame con una columna adicional 'status'
label_de_status <- function(results_df, lfc_threshold = 1, padj_threshold = 0.05) {
  results_df$status <- "NS"
  up <- !is.na(results_df$padj) &
    results_df$padj < padj_threshold &
    results_df$log2FoldChange >= lfc_threshold
  down <- !is.na(results_df$padj) &
    results_df$padj < padj_threshold &
    results_df$log2FoldChange <= -lfc_threshold

  results_df$status[up] <- "UP"
  results_df$status[down] <- "DOWN"
  results_df
}

#' Parseo simple de argumentos --clave valor desde la línea de comandos
#'
#' @param args vector de caracteres, típicamente commandArgs(trailingOnly = TRUE)
#' @param defaults lista nombrada con valores por defecto
#' @return lista nombrada con los argumentos resueltos
parse_cli_args <- function(args, defaults) {
  parsed <- defaults
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (key %in% names(defaults)) {
      parsed[[key]] <- args[i + 1]
      i <- i + 2
    } else {
      warning(sprintf("Argumento desconocido ignorado: %s", args[i]))
      i <- i + 1
    }
  }
  parsed
}
