#!/usr/bin/env Rscript
# run_pipeline.R
#
# Punto de entrada del pipeline de expresión diferencial.
# Uso:
#   Rscript run_pipeline.R --counts data/counts_raw.csv \
#                           --metadata data/sample_metadata.csv \
#                           --condition_col condition \
#                           --reference_level control \
#                           --outdir results/

suppressPackageStartupMessages({
  library(DESeq2)
})

script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) "."
)
if (is.null(script_dir) || script_dir == "") script_dir <- "."

source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "load_data.R"))
source(file.path(script_dir, "R", "qc_filtering.R"))
source(file.path(script_dir, "R", "deseq2_analysis.R"))
source(file.path(script_dir, "R", "visualization.R"))

defaults <- list(
  counts          = "data/counts_raw.csv",
  metadata        = "data/sample_metadata.csv",
  condition_col   = "condition",
  reference_level = "control",
  outdir          = "results",
  min_total_count = "10",
  lfc_threshold   = "1",
  padj_threshold  = "0.05"
)

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), defaults)

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

message(">> Cargando datos...")
data <- load_rnaseq_data(args$counts, args$metadata)

message(">> Filtrando genes de baja expresión...")
filtered_counts <- filter_low_expression(
  data$count_matrix,
  min_total_count = as.numeric(args$min_total_count)
)
message(sprintf(
  "   %d / %d genes retenidos tras el filtrado",
  nrow(filtered_counts), nrow(data$count_matrix)
))

message(">> Construyendo el modelo DESeq2...")
dds <- build_dds(
  filtered_counts, data$metadata,
  condition_col = args$condition_col,
  reference_level = args$reference_level
)

message(">> PCA exploratorio...")
run_exploratory_pca(
  dds, args$condition_col,
  out_path = file.path(args$outdir, "plots", "pca.png")
)

message(">> Ajustando DESeq2 y extrayendo resultados...")
results_df <- run_deseq2(dds, args$condition_col, args$reference_level)
results_df <- label_de_status(
  results_df,
  lfc_threshold = as.numeric(args$lfc_threshold),
  padj_threshold = as.numeric(args$padj_threshold)
)

write.csv(results_df, file.path(args$outdir, "deseq2_results.csv"), row.names = FALSE)
write.csv(
  subset(results_df, status != "NS"),
  file.path(args$outdir, "de_genes_significant.csv"),
  row.names = FALSE
)

message(">> Generando visualizaciones...")
plot_volcano(results_df, out_path = file.path(args$outdir, "plots", "volcano.png"))

# Igual que en qc_filtering.R: con pocos genes, vst() falla porque no
# alcanza el mínimo de submuestreo (nsub). Usamos la transformación
# completa, que funciona para cualquier tamaño de dataset.
vst_data <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
plot_top_variable_heatmap(
  SummarizedExperiment::assay(vst_data), data$metadata, args$condition_col,
  out_path = file.path(args$outdir, "plots", "heatmap_top_genes.png")
)

n_sig <- sum(results_df$status != "NS", na.rm = TRUE)
message(sprintf(">> Listo. %d genes significativos (padj < %s). Resultados en '%s'.",
                 n_sig, args$padj_threshold, args$outdir))
