# Small synthetic SpatialGraph fixture used across step-3 tests.
# Layout:
#   region "r1": cells a,b,c,d at the corners of a unit square
#   region "r2": cells e,f,g,h offset by (10,10), same square shape
# Each region is fully connected internally; one cross-region edge (d <-> e)
# lets us verify that subset_region drops cross-region edges.
fixture_small <- function() {
  nodes <- data.frame(
    cell_id   = c("a","b","c","d","e","f","g","h"),
    x         = c(0, 1, 0, 1, 10, 11, 10, 11),
    y         = c(0, 0, 1, 1, 10, 10, 11, 11),
    cell_type = c("T","B","T","B","T","B","T","B"),
    sample_id = c(rep("s1", 4), rep("s2", 4)),
    region_id = c(rep("r1", 4), rep("r2", 4)),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    source   = c("a","a","a","b","b","c","e","e","e","f","f","g","d"),
    target   = c("b","c","d","c","d","d","f","g","h","g","h","h","e"),
    distance = c(1, 1, sqrt(2), sqrt(2), 1, 1,
                 1, 1, sqrt(2), sqrt(2), 1, 1, sqrt(200)),
    sample_id = c(rep("s1", 6), rep("s2", 6), "s1"),
    stringsAsFactors = FALSE
  )
  list(nodes = nodes, edges = edges)
}

# Two-region fixture: spatially separated regions with different
# cell-type compositions so region-stratified shuffles preserve composition
# differently than global shuffles. Used by enrichment + compare tests.
fixture_two_regions <- function(N_per = 30L, seed = 42L) {
  set.seed(seed)
  r1 <- data.frame(
    cell_id   = paste0("r1_", seq_len(N_per)),
    x = stats::runif(N_per, 0, 1), y = stats::runif(N_per, 0, 1),
    cell_type = sample(c("T", "B"), N_per, replace = TRUE,
                       prob = c(0.85, 0.15)),
    sample_id = "s1", region_id = "r1",
    stringsAsFactors = FALSE
  )
  r2 <- data.frame(
    cell_id   = paste0("r2_", seq_len(N_per)),
    x = stats::runif(N_per, 5, 6), y = stats::runif(N_per, 5, 6),
    cell_type = sample(c("T", "B"), N_per, replace = TRUE,
                       prob = c(0.15, 0.85)),
    sample_id = "s1", region_id = "r2",
    stringsAsFactors = FALSE
  )
  rbind(r1, r2)
}

# Three-sample synthetic SpatialGraph + expression matrix for per-sample
# meta-analysis tests. Sample-level structure:
#   - 60 cells per sample, randomly typed T or B.
#   - Per-sample kNN graph (k = 4); no cross-sample edges.
#   - Expression: 30 genes, Poisson background.
#     * `g1` is enriched in T cells in ALL THREE samples — the
#       reproducible true positive.
#     * `g2` is enriched in T cells ONLY in sample s1 with a much
#       larger effect size — the "one-sample driver" pattern that
#       per-sample meta-analysis is supposed to flag via sd_log2FC.
fixture_multi_sample <- function(seed = 1L) {
  build_one <- function(sid, n, sample_seed) {
    set.seed(sample_seed)
    coords <- cbind(x = stats::runif(n), y = stats::runif(n))
    types <- sample(c("T", "B"), n, replace = TRUE)
    data.frame(
      cell_id   = paste0(sid, "_c", seq_len(n)),
      x         = coords[, "x"], y = coords[, "y"],
      cell_type = types,
      sample_id = sid,
      stringsAsFactors = FALSE
    )
  }
  s1 <- build_one("s1", 60L, seed)
  s2 <- build_one("s2", 60L, seed + 1L)
  s3 <- build_one("s3", 60L, seed + 2L)

  build_edges <- function(nodes_df) {
    sg <- build_spatial_graph(
      coords     = nodes_df[, c("x", "y")],
      cell_types = nodes_df$cell_type,
      sample_id  = nodes_df$sample_id[1L],
      cell_id    = nodes_df$cell_id,
      method     = "knn", k = 4L
    )
    sg$edges
  }
  e1 <- build_edges(s1); e2 <- build_edges(s2); e3 <- build_edges(s3)
  all_nodes <- rbind(s1, s2, s3)
  all_edges <- rbind(e1, e2, e3)
  sg <- SpatialGraph(all_nodes, all_edges,
                     meta = list(build_method = "synthetic_multi_sample"))

  set.seed(seed + 100L)
  G <- 30L; N <- nrow(all_nodes)
  expr <- matrix(stats::rpois(N * G, 0.5), nrow = G)
  rownames(expr) <- paste0("g", seq_len(G))
  colnames(expr) <- all_nodes$cell_id
  is_T <- all_nodes$cell_type == "T"

  # g1: reproducible across all three samples
  expr["g1", is_T] <- expr["g1", is_T] +
    stats::rpois(sum(is_T), 5)
  # g2: only enriched in s1, but with a larger spike
  in_s1_T <- all_nodes$sample_id == "s1" & is_T
  expr["g2", in_s1_T] <- expr["g2", in_s1_T] +
    stats::rpois(sum(in_s1_T), 12)

  list(sg = sg, expr = expr)
}

# Random scattered cells in the unit square, useful for kNN/Delaunay/radius
# tests. Seeded so failures are reproducible.
fixture_random <- function(N = 50L, seed = 1L) {
  set.seed(seed)
  data.frame(
    cell_id   = paste0("c", seq_len(N)),
    x         = runif(N),
    y         = runif(N),
    cell_type = sample(c("T","B","M"), N, replace = TRUE),
    sample_id = "s1",
    stringsAsFactors = FALSE
  )
}
