# rnaseq-de-explorer

Pipeline de análisis de expresión diferencial (RNA-seq) en R, basado en
DESeq2 y tidyverse, pensado para partir de una matriz de conteos y
llegar a una lista de genes diferencialmente expresados con visualizaciones
listas para interpretar (PCA, volcano plot, heatmap de los top genes).

No incluye el paso de cuantificación (alineamiento/pseudoalineamiento) —
asume que partimos una matriz de conteos por gen (p. ej. generada con
salmon, STAR+featureCounts, o kallisto). El foco de este pipeline es el análisis
estadístico y la interpretación de resultados.

## Estructura
```
.
├── run_pipeline.R          # Punto de entrada: ejecuta todo el análisis
├── R/
│   ├── load_data.R         # Carga y valida counts + metadata
│   ├── qc_filtering.R      # Filtrado de genes de baja expresión, PCA exploratorio
│   ├── deseq2_analysis.R   # Modelo DESeq2, contraste y tabla de resultados
│   ├── visualization.R     # Volcano plot, heatmap, MA plot
│   └── utils.R             # Funciones auxiliares reutilizadas por los scripts
├── data/
│   ├── counts_raw.csv      # Matriz de conteos de ejemplo (sintética)
│   └── sample_metadata.csv # Metadata de las muestras (condición, batch)
├── results/                # Salidas generadas por el pipeline (gitignored salvo .gitkeep)
├── tests/testthat/         # Tests unitarios de las funciones en R/utils.R
└── .github/workflows/ci.yml
```

## Requisitos
- R >= 4.2
- Paquetes: `DESeq2`, `tidyverse`, `pheatmap`, `ggrepel`, `testthat`

Instalación rápida de dependencias:
```r
install.packages(c("tidyverse", "pheatmap", "ggrepel", "testthat", "BiocManager"))
BiocManager::install("DESeq2")
```

## Uso:
```bash
Rscript run_pipeline.R \
  --counts data/counts_raw.csv \
  --metadata data/sample_metadata.csv \
  --condition_col condition \
  --reference_level control \
  --outdir results/
```

O interactivamente desde RStudio, abriendo `run_pipeline.R` y ejecutando por
bloques (útil para inspeccionar cada paso).

## Datos de ejemplo: 
`data/counts_raw.csv` 

## Salidas:
- `results/deseq2_results.csv`: tabla completa de resultados (log2FC, p-valor, padj) por gen.
- `results/de_genes_padj0.05.csv`: subconjunto de genes significativos (padj < 0.05).
- `results/plots/pca.png`: PCA de las muestras coloreado por condición.
- `results/plots/volcano.png`: volcano plot con genes significativos resaltados.
- `results/plots/heatmap_top_genes.png`: heatmap de los 20 genes más variables.
