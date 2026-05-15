# A SpatialGraph + expression fixture where one cell type ("T") strongly
# co-expresses a known gene, so a T-T-T motif's "in" cells should pop out
# on a target gene under the Wilcoxon test.
fixture_expr <- function(N = 60L, n_genes = 30L, seed = 11L) {
  set.seed(seed)
  fx <- fixture_random(N = N, seed = seed)
  sg <- build_spatial_graph(
    coords = fx[, c("x", "y")], cell_types = fx$cell_type,
    sample_id = fx$sample_id, cell_id = fx$cell_id,
    method = "knn", k = 4L
  )
  # Pure-Poisson background, then add a strong hit for "T" cells in g1.
  expr <- matrix(stats::rpois(N * n_genes, lambda = 0.5),
                 nrow = n_genes, ncol = N)
  rownames(expr) <- paste0("g", seq_len(n_genes))
  colnames(expr) <- sg$nodes$cell_id
  is_T <- sg$nodes$cell_type == "T"
  expr["g1", is_T] <- expr["g1", is_T] +
    stats::rpois(sum(is_T), lambda = 5)
  list(sg = sg, expr = expr,
       expr_sparse = Matrix::Matrix(expr, sparse = TRUE))
}

test_that("associate_genes(wilcoxon): output schema, finite stats, BH adjustment", {
  fx <- fixture_expr()
  m <- find_motifs(fx$sg, size = 3L)
  res <- associate_genes(m, fx$expr, fx$sg, method = "wilcoxon",
                         min_cells = 8L, min_in_expressed = 1L)
  expect_s3_class(res, "data.table")
  expect_setequal(
    colnames(res),
    c("motif_id", "orbit_id", "gene", "n_in", "n_out", "n_in_expressed",
      "mean_in", "mean_out", "log2FC", "U", "z", "p", "p_adj_BH")
  )
  expect_true(all(is.finite(res$z)))
  expect_true(all(res$p >= 0 & res$p <= 1))
  expect_true(all(res$p_adj_BH >= res$p))
})

test_that("the top-count motif gets at least one tested gene", {
  fx <- fixture_expr()
  m <- find_motifs(fx$sg, size = 3L)
  top <- m$catalog[which.max(m$catalog$count), motif_id]
  res <- associate_genes(m, fx$expr, fx$sg, method = "wilcoxon",
                         min_cells = 5L, min_in_expressed = 1L)
  expect_gte(sum(res$motif_id == top), 1L)
})

test_that("synthetic spike-in: g1 enriched in T-T-T motifs (z > 0, ranks high)", {
  fx <- fixture_expr(N = 80L, seed = 5L)
  m <- find_motifs(fx$sg, size = 3L)
  res <- associate_genes(m, fx$expr, fx$sg, method = "wilcoxon",
                         min_cells = 5L, min_in_expressed = 1L)
  # Find a T-only motif among the catalog (closed or open) with the most
  # cells, then check that g1 is among the lowest-p genes for it.
  ttt <- res[grepl("_T-T-T$|_T-T_T$", motif_id)]
  if (!nrow(ttt)) skip("no T-only size-3 motif present in this draw")
  best_t <- ttt[which.min(p)]
  expect_equal(best_t$gene, "g1")
  expect_gt(best_t$z, 0)
})

test_that("min_cells filter excludes small motifs", {
  fx <- fixture_expr(N = 30L, seed = 21L)
  m <- find_motifs(fx$sg, size = 3L)
  big_min <- 100L  # higher than any motif's cell count
  res <- associate_genes(m, fx$expr, fx$sg, method = "wilcoxon",
                         min_cells = big_min)
  expect_equal(nrow(res), 0L)
})

test_that("by = 'orbit' yields one row per (motif, orbit_id, gene)", {
  fx <- fixture_expr(N = 60L, seed = 3L)
  m <- find_motifs(fx$sg, size = 3L)
  res_motif <- associate_genes(m, fx$expr, fx$sg, by = "motif",
                                min_cells = 6L, min_in_expressed = 1L)
  res_orbit <- associate_genes(m, fx$expr, fx$sg, by = "orbit",
                                min_cells = 6L, min_in_expressed = 1L)
  expect_true(all(is.na(res_motif$orbit_id)))
  expect_false(all(is.na(res_orbit$orbit_id)))
  # Orbit-level testing should produce at least as many rows as motif-level
  # for any motif that has > 1 orbit.
  by_motif <- res_motif[, .(n = .N), by = motif_id]
  by_orbit <- res_orbit[, .(n = .N), by = motif_id]
  joined <- merge(by_motif, by_orbit, by = "motif_id",
                  suffixes = c("_motif", "_orbit"))
  expect_true(all(joined$n_orbit >= joined$n_motif))
})

test_that("module_score: schema check + module names round-trip into 'gene' column", {
  fx <- fixture_expr(N = 60L, seed = 8L)
  m <- find_motifs(fx$sg, size = 3L)
  modules <- list(
    spike   = c("g1"),
    random5 = paste0("g", 2:6)
  )
  res <- associate_genes(m, fx$expr, fx$sg, method = "module_score",
                         modules = modules, min_cells = 6L)
  expect_s3_class(res, "data.table")
  expect_true(all(res$gene %in% c("spike", "random5")))
  expect_true(all(is.finite(res$z)))
})

test_that("sparse and dense expression give identical Wilcoxon results", {
  fx <- fixture_expr(N = 50L, seed = 17L)
  m <- find_motifs(fx$sg, size = 3L)
  r_dense  <- associate_genes(m, fx$expr,        fx$sg, min_cells = 5L,
                               min_in_expressed = 1L)
  r_sparse <- associate_genes(m, fx$expr_sparse, fx$sg, min_cells = 5L,
                               min_in_expressed = 1L)
  data.table::setkey(r_dense,  motif_id, orbit_id, gene)
  data.table::setkey(r_sparse, motif_id, orbit_id, gene)
  expect_equal(r_dense$U, r_sparse$U)
  expect_equal(r_dense$p, r_sparse$p)
})

test_that("per_sample mode produces the meta-analysis schema", {
  fx <- fixture_multi_sample()
  m <- find_motifs(fx$sg, size = 3L)
  res <- associate_genes(
    m, fx$expr, fx$sg, method = "wilcoxon",
    min_cells = 5L, min_in_expressed = 1L,
    pooling = "per_sample"
  )
  expect_s3_class(res, "data.table")
  expect_setequal(
    colnames(res),
    c("motif_id", "orbit_id", "gene", "n_samples_tested",
      "n_in_total", "n_out_total", "mean_log2FC", "sd_log2FC",
      "z_combined", "p_combined", "p_adj_BH")
  )
  expect_true(all(res$p_combined >= 0 & res$p_combined <= 1))
  expect_true(all(is.finite(res$z_combined)))
  expect_true(all(res$n_samples_tested >= 1L))
  expect_true(all(res$p_adj_BH >= res$p_combined))
})

test_that("reproducible gene g1 has lower sd_log2FC than heterogeneous g2", {
  fx <- fixture_multi_sample()
  m <- find_motifs(fx$sg, size = 3L)
  res <- associate_genes(
    m, fx$expr, fx$sg, method = "wilcoxon",
    min_cells = 5L, min_in_expressed = 1L,
    pooling = "per_sample"
  )
  # Restrict to T-only motifs where the spike-ins live (largest "in" set).
  ttt <- res[grepl("_T-T-T$|_T-T_T$", motif_id) &
             n_samples_tested >= 2L]
  if (!nrow(ttt[gene == "g1"]) || !nrow(ttt[gene == "g2"])) {
    skip("not enough T-only motifs sampled across multiple samples")
  }
  best_g1 <- ttt[gene == "g1"][which.min(p_combined)]
  best_g2 <- ttt[gene == "g2"][which.min(p_combined)]
  # Both should be nominally enriched (positive z_combined).
  expect_gt(best_g1$z_combined, 0)
  expect_gt(best_g2$z_combined, 0)
  # The reproducible gene should have a smaller spread of log2FCs.
  expect_lt(best_g1$sd_log2FC, best_g2$sd_log2FC)
})

test_that("per_sample on a single-sample sg falls back to pooled with a warning", {
  fx <- fixture_random(N = 30L)
  sg <- build_spatial_graph(
    fx[, c("x", "y")], fx$cell_type, fx$sample_id,
    cell_id = fx$cell_id, method = "knn", k = 3L
  )
  expr <- matrix(stats::rpois(30L * 10L, 1), nrow = 10L)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- sg$nodes$cell_id
  m <- find_motifs(sg, size = 3L)
  expect_warning(
    res <- associate_genes(m, expr, sg, method = "wilcoxon",
                           pooling = "per_sample",
                           min_cells = 5L, min_in_expressed = 1L),
    "falls back to 'pooled'"
  )
  # Schema should match the pooled output, not the per-sample output.
  expect_true("z" %in% colnames(res))
  expect_true("p" %in% colnames(res))
  expect_false("z_combined" %in% colnames(res))
})

test_that("argument validation: errors on bad inputs", {
  fx <- fixture_expr(N = 20L, seed = 4L)
  m <- find_motifs(fx$sg, size = 3L)
  expect_error(associate_genes(m, "not a matrix", fx$sg), "matrix")
  bad <- fx$expr; rownames(bad) <- NULL
  expect_error(associate_genes(m, bad, fx$sg), "row names")
  expect_error(
    associate_genes(m, fx$expr, fx$sg, method = "module_score"),
    "modules"
  )
})
