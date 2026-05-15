test_that("constructor validates required columns", {
  expect_error(
    SpatialGraph(data.frame(x = 1), data.frame(source = "a", target = "b",
                                                sample_id = "s")),
    "nodes is missing required columns"
  )
  expect_error(
    SpatialGraph(
      data.frame(cell_id = "a", x = 0, y = 0, cell_type = "T",
                 sample_id = "s", stringsAsFactors = FALSE),
      data.frame(source = "a", target = "b")
    ),
    "edges is missing required columns"
  )
})

test_that("constructor canonicalizes edges (source < target, dedup, no self-loops)", {
  fx <- fixture_small()
  # Inject a swapped edge and a duplicate to verify normalization.
  bad_edges <- rbind(
    fx$edges,
    data.frame(source = "b", target = "a", distance = 1,
               sample_id = "s1", stringsAsFactors = FALSE),
    data.frame(source = "a", target = "a", distance = 0,
               sample_id = "s1", stringsAsFactors = FALSE)
  )
  expect_warning(
    sg <- SpatialGraph(fx$nodes, bad_edges),
    "self-loop"
  )
  expect_true(all(sg$edges$source < sg$edges$target))
  expect_equal(anyDuplicated(sg$edges, by = c("source","target")), 0L)
})

test_that("constructor errors on duplicate cell_ids and dangling endpoints", {
  fx <- fixture_small()
  dup <- rbind(fx$nodes, fx$nodes[1, ])
  expect_error(SpatialGraph(dup, fx$edges), "duplicates")

  bad_edges <- rbind(
    fx$edges,
    data.frame(source = "a", target = "z_not_real",
               distance = 1, sample_id = "s1",
               stringsAsFactors = FALSE)
  )
  expect_error(SpatialGraph(fx$nodes, bad_edges), "endpoint")
})

test_that("read -> write -> read produces an identical SpatialGraph", {
  fx <- fixture_small()
  sg1 <- SpatialGraph(fx$nodes, fx$edges, meta = list(build_method = "manual"))
  td <- withr::local_tempdir()
  np <- file.path(td, "nodes.parquet")
  ep <- file.path(td, "edges.parquet")
  write_spatial_graph(sg1, np, ep)
  sg2 <- read_spatial_graph(np, ep)
  np2 <- file.path(td, "nodes2.parquet")
  ep2 <- file.path(td, "edges2.parquet")
  write_spatial_graph(sg2, np2, ep2)
  sg3 <- read_spatial_graph(np2, ep2)

  expect_equal(sg2$nodes, sg3$nodes)
  expect_equal(sg2$edges, sg3$edges)
  # File-path keys differ across the two reads; only structural meta should match.
  expect_equal(sg2$meta$build_method, sg3$meta$build_method)
})

test_that("subset_region keeps only edges with both endpoints inside", {
  fx <- fixture_small()
  sg <- SpatialGraph(fx$nodes, fx$edges)
  sg_r1 <- subset_region(sg, region_id == "r1")
  expect_lt(nrow(sg_r1$nodes), nrow(sg$nodes))
  ids <- sg_r1$nodes$cell_id
  expect_true(all(sg_r1$edges$source %in% ids))
  expect_true(all(sg_r1$edges$target %in% ids))
  # The cross-region edge d <-> e must be gone.
  cross <- sg_r1$edges[(source == "d" & target == "e") |
                       (source == "e" & target == "d")]
  expect_equal(nrow(cross), 0L)
})

test_that("subset_sample composes with subset_region", {
  fx <- fixture_small()
  sg <- SpatialGraph(fx$nodes, fx$edges)
  sg2 <- subset_sample(sg, sample_id == "s1")
  expect_true(all(sg2$nodes$sample_id == "s1"))
  expect_true(all(sg2$edges$source %in% sg2$nodes$cell_id))
})

test_that("as_igraph builds the right vertex and edge counts and caches the result", {
  fx <- fixture_small()
  sg <- SpatialGraph(fx$nodes, fx$edges)
  g  <- as_igraph(sg)
  expect_s3_class(g, "igraph")
  expect_equal(igraph::vcount(g), nrow(sg$nodes))
  expect_equal(igraph::ecount(g), nrow(sg$edges))
  expect_false(igraph::is_directed(g))
  expect_identical(as_igraph(sg), g)  # same object; cache hit
})

test_that("as_igraph cache invalidates when the underlying SpatialGraph changes", {
  fx <- fixture_small()
  sg <- SpatialGraph(fx$nodes, fx$edges)
  g1 <- as_igraph(sg)
  sg2 <- subset_region(sg, region_id == "r1")
  g2 <- as_igraph(sg2)
  expect_lt(igraph::vcount(g2), igraph::vcount(g1))
})
