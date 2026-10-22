# smotif 0.2.1

- Every error for a missing smotifrs now gives the R-universe install line
  (pre-built binaries, no Rust toolchain). Only `find_motifs(backend =
  "rust")` did before; `motif_enrichment_stream()`,
  `motif_enrichment_store()`, `motif_instances()` and `store_node_order()`
  only said smotifrs was needed. When smotifrs is installed but
  `backend = "r"` was chosen, the size-4 and conditional-null refusals point
  at `backend = "auto"` or `"rust"` instead.

# smotif 0.2.0

- `motif_enrichment_stream()`: enrichment over an Arrow stream of edges, with
  labels passed as a lookup keyed by integer node id. The network's nodes are
  the stream's endpoints. Needs smotifrs >= 0.2.0.

# smotif 0.1.0 (in development)

Initial scaffold of the `smotif` package for spatial multi-cellular motif
discovery in spatial omics data.

## In scope for v0.1

- `SpatialGraph` S3 class with parquet round-trip via `arrow`.
- `build_spatial_graph()` over coordinates with kNN, Delaunay, and radius methods.
- Motif enumeration (size 2 and 3 closed/open, size 4 via `igraph::motifs`)
  with colored canonical labeling for cell-typed nodes.
- Spatially-aware null models: label permutation, region-stratified label
  permutation, geometric jitter.
- Empirical enrichment testing with BH adjustment.
- Wilcoxon and module-score gene association at motif and orbit level.
- Cross-region and cross-sample motif comparison.
- Reference dataset preparation script (`inst/scripts/prepare_workshop_xenium.R`)
  for the workshop Xenium human lung cancer FFPE slide.

## Out of scope for v0.1 (planned later)

- Rust / extendr backend (a `# TODO(rust-backend)` marker in `R/motifs.R`
  flags the dispatch site).
- Motifs of size 5 or larger.
- Directed motifs, hypergraphs, motif embeddings.
- Per-sample meta-analysis in `associate_genes()`.
- Visualization beyond `plot.Motifs()`.
- Giotto Suite integration (intended as a separate downstream package).
- Full geometric resampling null with re-thresholding.
