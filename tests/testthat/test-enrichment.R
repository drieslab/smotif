test_that("test_motif_enrichment output schema is well-formed", {
  fx <- fixture_random(N = 60L)
  sg <- build_spatial_graph(fx[, c("x", "y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 4L)
  m <- find_motifs(sg, size = 3L)
  res <- test_motif_enrichment(m, sg, null = "label_permutation",
                               n_perm = 100L, seed = 1L)
  expect_s3_class(res, "data.table")
  expect_setequal(
    colnames(res),
    c("motif_id", "observed", "mean_null", "sd_null", "z",
      "p_emp", "p_adj_BH", "fold")
  )
  expect_equal(nrow(res), nrow(m$catalog))
  expect_true(all(res$p_emp > 0 & res$p_emp <= 1))
  expect_true(all(is.finite(res$z)))
  expect_true(all(res$p_adj_BH >= res$p_emp))
})

test_that("region_stratified differs from label_permutation when regions differ", {
  fx <- fixture_two_regions(N_per = 40L)
  sg <- build_spatial_graph(
    fx[, c("x", "y")], fx$cell_type, fx$sample_id,
    cell_id = fx$cell_id, method = "knn", k = 4L,
    region_id = fx$region_id
  )
  m <- find_motifs(sg, size = 3L)
  r_global <- test_motif_enrichment(m, sg, null = "label_permutation",
                                     n_perm = 200L, seed = 1L)
  r_strat  <- test_motif_enrichment(m, sg, null = "region_stratified",
                                     n_perm = 200L, seed = 1L)
  data.table::setkey(r_global, motif_id)
  data.table::setkey(r_strat, motif_id)
  joined <- merge(r_global[, .(motif_id, p_global = p_emp,
                                m_global = mean_null)],
                  r_strat[, .(motif_id, p_strat = p_emp,
                               m_strat = mean_null)],
                  by = "motif_id")
  # Stratification should change the null mean for at least some motifs
  # (composition-preserving vs composition-mixing).
  expect_true(any(abs(joined$m_global - joined$m_strat) > 0.5))
})

test_that("region_stratified errors without region_id", {
  fx <- fixture_random(N = 20L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 3L)
  # Strip region_id (build_spatial_graph leaves it absent unless passed).
  expect_false("region_id" %in% names(sg$nodes))
  m <- find_motifs(sg, size = 3L)
  expect_error(
    test_motif_enrichment(m, sg, null = "region_stratified", n_perm = 10L),
    "region_id"
  )
})

test_that("p_emp is reproducible under the same seed", {
  fx <- fixture_random(N = 30L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 4L)
  m <- find_motifs(sg, size = 3L)
  r1 <- test_motif_enrichment(m, sg, n_perm = 50L, seed = 7L)
  r2 <- test_motif_enrichment(m, sg, n_perm = 50L, seed = 7L)
  expect_equal(r1$p_emp, r2$p_emp)
  expect_equal(r1$mean_null, r2$mean_null)
})

test_that("geometric_jitter null runs and produces valid statistics", {
  skip_if_not_installed("FNN")
  fx <- fixture_random(N = 30L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 4L)
  m <- find_motifs(sg, size = 3L)
  res <- test_motif_enrichment(
    m, sg, null = "geometric_jitter", n_perm = 25L,
    seed = 1L, null_args = list(sd = 0.05)
  )
  expect_true(all(res$p_emp > 0 & res$p_emp <= 1))
  expect_true(all(is.finite(res$z)))
})

test_that("geometric_jitter requires a build_method we can rebuild from", {
  # Hand-built SpatialGraph carries no build_method we can use, so the
  # jitter null should refuse rather than silently produce something
  # nonsensical.
  fx <- fixture_small()
  sg <- SpatialGraph(fx$nodes, fx$edges, meta = list(build_method = "manual"))
  m <- find_motifs(sg, size = 2L)
  expect_error(
    test_motif_enrichment(m, sg, null = "geometric_jitter",
                          n_perm = 5L, null_args = list(sd = 0.01)),
    "build_method"
  )
})

test_that("uncolored permutation null preserves total count invariant", {
  fx <- fixture_random(N = 40L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 4L)
  m <- find_motifs(sg, size = 2L, colored = FALSE)
  res <- test_motif_enrichment(m, sg, n_perm = 50L, seed = 2L)
  # Without coloring, label permutation cannot change the count of
  # "size2_edge". Mean and observed must coincide exactly.
  expect_equal(res$mean_null, res$observed)
  expect_equal(res$sd_null, 0)
})
