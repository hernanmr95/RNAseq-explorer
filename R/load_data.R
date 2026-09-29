# R/load_data.R
# Carga la matriz de counts y la metadata de muestras, y valida que sean
# consistentes entre sí antes de pasar al análisis.

suppressPackageStartupMessages(library(dplyr))

#' Carga counts + metadata desde disco
#'
#' @param counts_path ruta al CSV de counts (gene_id en la primera columna)
#' @param metadata_path ruta al CSV de metadata (debe incluir sample_id)
#' @return lista con: $count_matrix (matriz numérica gene x sample) y $metadata (data.frame)
load_rnaseq_data <- function(counts_path, metadata_path) {
  if (!file.exists(counts_path)) {
    stop(sprintf("No se encuentra el archivo de counts: %s", counts_path))
  }
  if (!file.exists(metadata_path)) {
    stop(sprintf("No se encuentra el archivo de metadata: %s", metadata_path))
  }

  counts_df <- read.csv(counts_path, stringsAsFactors = FALSE, check.names = FALSE)
  metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

  if (!"sample_id" %in% colnames(metadata)) {
    stop("La metadata debe incluir una columna 'sample_id'")
  }

  validate_inputs(counts_df, metadata)

  count_matrix <- as.matrix(counts_df[, -1])
  rownames(count_matrix) <- counts_df[[1]]
  mode(count_matrix) <- "integer"

  # Reordenar metadata para que coincida exactamente con las columnas de counts
  metadata <- metadata[match(colnames(count_matrix), metadata$sample_id), ]

  list(count_matrix = count_matrix, metadata = metadata)
}
