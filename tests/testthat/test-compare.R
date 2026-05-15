test_that("compare_motifs returns one row per motif present in any group", {
  fx <- fixture_two_regions(N_per = 30L)
  sg <- build_spatial_graph(
    fx[, c("x", "y")], fx$cell_type, fx$sample_id,
    cell_id = fx$cell_id, region_id = fx$region_id,
    method = "knn", k = 4L
  )
  sa <- subset_region(sg, region_id == "r1")
  sb <- subset_region(sg, region_id == "r2")
  ma <- find_motifs(sa, size = 3L)
  mb <- find_motifs(sb, size = 3L)

  res <- compare_motifs(list(r1 = ma, r2 = mb))
  expect_s3_class(res, "data.table")
  expect_setequal(
    intersect(colnames(res),
              c("motif_id", "size", "canonical_iso_class", "color_tuple",
                "count_r1", "count_r2", "total", "chi2", "df",
                "p", "p_adj_BH")),
    c("motif_id", "size", "canonical_iso_class", "color_tuple",
      "count_r1", "count_r2", "total", "chi2", "df", "p", "p_adj_BH")
  )

  union_motifs <- union(ma$catalog$motif_id, mb$catalog$motif_id)
  expect_setequal(res$motif_id, union_motifs)
  expect_true(all(res$total ==
                  res$count_r1 + res$count_r2))
  expect_true(all(res$p >= 0 & res$p <= 1))
  expect_true(all(res$p_adj_BH >= res$p))
})

test_that("motifs absent from one group keep zero counts there", {
  fx <- fixture_two_regions(N_per = 30L)
  sg <- build_spatial_graph(
    fx[, c("x", "y")], fx$cell_type, fx$sample_id,
    cell_id = fx$cell_id, region_id = fx$region_id,
    method = "knn", k = 4L
  )
  sa <- subset_region(sg, region_id == "r1")
  sb <- subset_region(sg, region_id == "r2")
  ma <- find_motifs(sa, size = 3L)
  mb <- find_motifs(sb, size = 3L)
  res <- compare_motifs(list(r1 = ma, r2 = mb))

  only_r1 <- setdiff(ma$catalog$motif_id, mb$catalog$motif_id)
  if (length(only_r1)) {
    sub <- res[motif_id %in% only_r1]
    expect_true(all(sub$count_r2 == 0L))
    expect_true(all(sub$count_r1 > 0L))
  }
})

test_that("3-niche comparison runs and df = n_groups - 1", {
  set.seed(2)
  N <- 90L
  nodes <- data.frame(
    cell_id   = paste0("c", seq_len(N)),
    x         = stats::runif(N), y = stats::runif(N),
    cell_type = sample(c("T", "B"), N, replace = TRUE),
    sample_id = "s1",
    region_id = sample(c("a", "b", "c"), N, replace = TRUE),
    stringsAsFactors = FALSE
  )
  sg <- build_spatial_graph(nodes[, c("x","y")], nodes$cell_type,
                            nodes$sample_id, cell_id = nodes$cell_id,
                            region_id = nodes$region_id,
                            method = "knn", k = 4L)
  ml <- list(
    a = find_motifs(subset_region(sg, region_id == "a"), size = 3L),
    b = find_motifs(subset_region(sg, region_id == "b"), size = 3L),
    c = find_motifs(subset_region(sg, region_id == "c"), size = 3L)
  )
  res <- compare_motifs(ml)
  expect_true(all(res$df == 2L))
  expect_true(all(c("count_a", "count_b", "count_c") %in% colnames(res)))
})

test_that("group argument lets replicates collapse into named groups", {
  fx <- fixture_two_regions(N_per = 25L)
  sg <- build_spatial_graph(
    fx[, c("x","y")], fx$cell_type, fx$sample_id,
    cell_id = fx$cell_id, region_id = fx$region_id,
    method = "knn", k = 3L
  )
  sa <- subset_region(sg, region_id == "r1")
  sb <- subset_region(sg, region_id == "r2")
  ma1 <- find_motifs(sa, size = 2L)
  ma2 <- find_motifs(sa, size = 2L)
  mb  <- find_motifs(sb, size = 2L)
  res <- compare_motifs(
    list(a1 = ma1, a2 = ma2, b = mb),
    group = c("A", "A", "B")
  )
  expect_true(all(c("count_A", "count_B") %in% colnames(res)))
  expect_false(any(c("count_a1", "count_a2", "count_b") %in% colnames(res)))
  # Replicates should sum within their group.
  any_motif <- res$motif_id[1]
  combined <- ma1$catalog[motif_id == any_motif]$count +
              ma2$catalog[motif_id == any_motif]$count
  if (length(combined) == 0L) combined <- 0L
  expect_equal(res[motif_id == any_motif, count_A], combined)

  expect_error(compare_motifs(list(only = ma1)), "two groups")
})

test_that("plot.Motifs returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  fx <- fixture_random(N = 30L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 4L)
  m <- find_motifs(sg, size = 3L)
  p <- plot(m, top = 5L)
  expect_s3_class(p, "ggplot")
})

test_that("compare_motifs validates inputs", {
  expect_error(compare_motifs(list()), "non-empty")
  expect_error(compare_motifs(list(a = "not_motifs")), "Motifs")
})
