# tests/testthat/test-utils.R

source(file.path("..", "..", "R", "utils.R"))

test_that("filter_low_expression elimina genes por debajo del umbral", {
  mat <- matrix(
    c(5, 2, 1, 50, 60, 55, 0, 0, 1),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("geneA", "geneB", "geneC"), c("s1", "s2", "s3"))
  )
  filtered <- filter_low_expression(mat, min_total_count = 10)

  expect_equal(rownames(filtered), "geneB")
  expect_equal(nrow(filtered), 1)
})

test_that("label_de_status clasifica correctamente UP/DOWN/NS", {
  df <- data.frame(
    gene_id = c("g1", "g2", "g3", "g4"),
    log2FoldChange = c(2.5, -3.0, 0.2, 1.5),
    padj = c(0.001, 0.01, 0.02, 0.2)
  )
  labeled <- label_de_status(df, lfc_threshold = 1, padj_threshold = 0.05)

  expect_equal(labeled$status, c("UP", "DOWN", "NS", "NS"))
})

test_that("validate_inputs lanza error si hay muestras que no coinciden", {
  counts <- data.frame(gene_id = "g1", s1 = 10, s2 = 20)
  metadata <- data.frame(sample_id = c("s1", "s3"), condition = c("control", "treated"))

  expect_error(validate_inputs(counts, metadata), "no están en la metadata|no en la metadata")
})

test_that("parse_cli_args resuelve argumentos y mantiene defaults no provistos", {
  defaults <- list(a = "1", b = "2")
  parsed <- parse_cli_args(c("--a", "99"), defaults)

  expect_equal(parsed$a, "99")
  expect_equal(parsed$b, "2")
})
