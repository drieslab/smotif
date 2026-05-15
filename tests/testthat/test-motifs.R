# Brute-force enumerator used as ground truth for the canonical-labeling
# sanity check. Walks all 3-vertex subsets and bins them by topology and
# sorted color tuple. O(N^3); only used on tiny fixtures.
brute_size3 <- function(nodes, edges) {
  ids <- nodes$cell_id
  ct  <- stats::setNames(nodes$cell_type, nodes$cell_id)
  ekey <- paste(pmin(edges$source, edges$target),
                pmax(edges$source, edges$target), sep = "|")
  has_edge <- function(a, b) {
    paste(min(a, b), max(a, b), sep = "|") %in% ekey
  }
  out <- list()
  N <- length(ids)
  for (i in seq_len(N - 2L)) {
    for (j in seq.int(i + 1L, N - 1L)) {
      for (k in seq.int(j + 1L, N)) {
        a <- ids[i]; b <- ids[j]; c <- ids[k]
        ab <- has_edge(a, b); ac <- has_edge(a, c); bc <- has_edge(b, c)
        n_e <- ab + ac + bc
        if (n_e == 3L) {
          tag <- paste0("size3_closed_",
                        paste(sort(c(ct[a], ct[b], ct[c])), collapse = "-"))
          out[[length(out) + 1L]] <- tag
        } else if (n_e == 2L) {
          # Identify the center (the vertex with degree 2 in this subgraph).
          deg_a <- ab + ac; deg_b <- ab + bc; deg_c <- ac + bc
          centers <- c(a, b, c)[c(deg_a, deg_b, deg_c) == 2L]
          ends    <- c(a, b, c)[c(deg_a, deg_b, deg_c) == 1L]
          end_colors <- sort(ct[ends])
          tag <- paste0("size3_open_",
                        end_colors[1], "-", end_colors[2], "_", ct[centers])
          out[[length(out) + 1L]] <- tag
        }
      }
    }
  }
  tab <- table(unlist(out))
  dt <- data.table::data.table(
    motif_id = names(tab),
    count    = as.integer(tab)
  )
  data.table::setorderv(dt, "motif_id")
  dt
}

test_that("size 2: total count equals nrow(edges); colored partitions it", {
  fx <- fixture_random(N = 30L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 4L)
  m_un <- find_motifs(sg, size = 2L, colored = FALSE)
  expect_equal(nrow(m_un$catalog), 1L)
  expect_equal(sum(m_un$catalog$count), nrow(sg$edges))

  m <- find_motifs(sg, size = 2L, colored = TRUE)
  expect_equal(sum(m$catalog$count), nrow(sg$edges))
  # Color tuples are non-empty when colored = TRUE.
  expect_false(any(is.na(m$catalog$color_tuple)))
  # Every motif_id starts with "size2_edge_".
  expect_true(all(grepl("^size2_edge", m$catalog$motif_id)))
  # Two incidence rows per instance.
  expect_equal(nrow(m$incidence),
               2L * sum(m$catalog$count))
})

test_that("size 3 closed count matches igraph::triangles()", {
  fx <- fixture_random(N = 60L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 4L)
  g <- as_igraph(sg)
  n_tri <- length(igraph::triangles(g)) %/% 3L
  m <- find_motifs(sg, size = 3L, colored = FALSE)
  closed <- m$catalog[canonical_iso_class == "triangle"]
  expect_equal(sum(closed$count), n_tri)
})

test_that("size 3 motif catalog matches a brute-force enumerator (colored)", {
  fx <- fixture_random(N = 25L, seed = 7L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 3L)
  m <- find_motifs(sg, size = 3L, colored = TRUE)
  bf <- brute_size3(sg$nodes, sg$edges)

  agg <- m$catalog[, .(count = sum(count)), by = motif_id][order(motif_id)]
  expect_equal(agg$motif_id, bf$motif_id)
  expect_equal(agg$count, bf$count)
})

test_that("orbit_ids are consistent: vertices with same (deg,color) share orbit", {
  fx <- fixture_random(N = 40L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 3L)
  m <- find_motifs(sg, size = 3L, colored = TRUE)
  # Within a single instance, each orbit_id appears at most as many times
  # as there are matching (deg, color) vertices. Globally: for each
  # motif_id, all instances share an orbit-id histogram.
  inc <- merge(m$incidence, sg$nodes[, .(cell_id, cell_type)],
               by = "cell_id")
  per_inst <- inc[, .(sig = paste(sort(orbit_id), collapse = ",")),
                  by = .(motif_id, instance_id)]
  per_motif <- per_inst[, .(uniq = data.table::uniqueN(sig)),
                        by = motif_id]
  expect_true(all(per_motif$uniq == 1L))
})

test_that("anchored_on restricts enumeration to motifs touching at least one anchor", {
  fx <- fixture_random(N = 40L, seed = 3L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 4L)
  anchors <- fx$cell_id[1:3]
  m <- find_motifs(sg, size = 3L, anchored_on = anchors, colored = FALSE)
  # Every enumerated instance must include at least one anchored cell_id.
  by_inst <- m$incidence[, .(touches = any(cell_id %in% anchors)),
                         by = instance_id]
  expect_true(all(by_inst$touches))
  # Smoke-check that the unrestricted enumeration is at least as large.
  m_full <- find_motifs(sg, size = 3L, colored = FALSE)
  expect_gte(sum(m_full$catalog$count), sum(m$catalog$count))
})

test_that("max_instances caps total instances and preserves catalog/incidence consistency", {
  fx <- fixture_random(N = 60L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 4L)
  m_full <- find_motifs(sg, size = 3L, colored = FALSE)
  cap <- max(1L, sum(m_full$catalog$count) %/% 4L)
  expect_warning(
    m <- find_motifs(sg, size = 3L, colored = FALSE, max_instances = cap),
    "subsampling"
  )
  expect_lte(sum(m$catalog$count), cap)
  expect_equal(data.table::uniqueN(m$incidence$instance_id),
               sum(m$catalog$count))
})

test_that("rust backend dispatches to smotifrs when present, size validation always", {
  fx <- fixture_random(N = 10L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 3L)
  if (requireNamespace("smotifrs", quietly = TRUE)) {
    m <- find_motifs(sg, size = 3L, backend = "rust")
    expect_s3_class(m, "Motifs")
    expect_equal(m$meta$backend, "rust")
  } else {
    expect_error(find_motifs(sg, size = 3L, backend = "rust"),
                 "smotifrs")
  }
  expect_error(find_motifs(sg, size = 5L), "must be one of")
})

test_that("Motifs object prints without error", {
  fx <- fixture_random(N = 20L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                            cell_id = fx$cell_id, method = "knn", k = 3L)
  m <- find_motifs(sg, size = 3L)
  expect_output(print(m), "Motifs")
})
