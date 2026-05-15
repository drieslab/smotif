test_that("build_spatial_graph(method='knn') returns a valid SpatialGraph", {
  skip_if_not_installed("FNN")
  fx <- fixture_random(N = 50L)
  sg <- build_spatial_graph(
    coords = fx[, c("x", "y")],
    cell_types = fx$cell_type,
    sample_id = fx$sample_id,
    cell_id = fx$cell_id,
    method = "knn",
    k = 4L
  )
  expect_s3_class(sg, "SpatialGraph")
  expect_equal(nrow(sg$nodes), 50L)
  expect_true(all(sg$edges$source < sg$edges$target))
  # k=4 with undirected dedup gives at most k*N/2 ... 2*k*N/2 edges depending
  # on symmetry. Bound generously and check both ends.
  expect_gte(nrow(sg$edges), 50L * 4L / 2L)
  expect_lte(nrow(sg$edges), 50L * 4L)
  expect_equal(sg$meta$build_method, "knn")
  expect_equal(sg$meta$k, 4L)
})

test_that("build_spatial_graph respects max_distance", {
  skip_if_not_installed("FNN")
  fx <- fixture_random(N = 50L)
  sg_no_cap <- build_spatial_graph(
    fx[, c("x", "y")], fx$cell_type, "s1", method = "knn", k = 6L
  )
  cap <- median(sg_no_cap$edges$distance)
  sg_cap <- build_spatial_graph(
    fx[, c("x", "y")], fx$cell_type, "s1", method = "knn", k = 6L,
    max_distance = cap
  )
  expect_lt(nrow(sg_cap$edges), nrow(sg_no_cap$edges))
  expect_true(all(sg_cap$edges$distance <= cap))
})

test_that("build_spatial_graph(method='delaunay') runs and gives planar edge counts", {
  skip_if_not_installed("deldir")
  fx <- fixture_random(N = 30L)
  sg <- build_spatial_graph(
    fx[, c("x","y")], fx$cell_type, fx$sample_id, fx$cell_type,
    method = "delaunay"
  )
  # A planar triangulation on N points has at most 3N-6 edges.
  expect_lte(nrow(sg$edges), 3L * 30L - 6L)
  expect_equal(sg$meta$build_method, "delaunay")
})

test_that("build_spatial_graph(method='radius') matches a brute-force pairwise check", {
  skip_if_not_installed("dbscan")
  fx <- fixture_random(N = 25L)
  rad <- 0.25
  sg <- build_spatial_graph(
    fx[, c("x","y")], fx$cell_type, fx$sample_id,
    method = "radius", radius = rad
  )
  # Brute force on this small fixture.
  d <- as.matrix(dist(fx[, c("x", "y")]))
  pairs <- which(d > 0 & d <= rad, arr.ind = TRUE)
  pairs <- pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]
  expect_equal(nrow(sg$edges), nrow(pairs))
})

test_that("build_spatial_graph errors on missing required arguments", {
  fx <- fixture_random(N = 10L)
  expect_error(
    build_spatial_graph(fx[, c("x","y")], fx$cell_type[1:5], "s1",
                       method = "knn", k = 3L),
    "length\\(cell_types\\)"
  )
  expect_error(
    build_spatial_graph(fx[, c("x","y")], fx$cell_type, "s1",
                       method = "radius"),
    "requires a non-null `radius`"
  )
})
