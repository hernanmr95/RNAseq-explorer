# Changelog

### Planeado
- Anotación funcional de los genes significativos (clusterProfiler / enrichR).
- Soporte para diseños con más de dos condiciones (contrastes múltiples).
- Corrección de batch en el diseño del modelo (`~ batch + condition`).

## [1.0.0] - 2026-07-10

### Añadido
- Pipeline completo: carga → filtrado → DESeq2 → visualización.
- Dataset sintético de ejemplo (40 genes x 6 muestras) para validar el flujo end-to-end.
- Tests unitarios de las funciones auxiliares (`tests/testthat/`).
