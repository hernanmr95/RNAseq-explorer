#!/usr/bin/env Rscript
# tests/run_tests.R
# Ejecuta todos los tests unitarios del proyecto.

suppressPackageStartupMessages(library(testthat))

test_dir("tests/testthat", reporter = "summary")
