# Brute-force enumerator for connected 4-vertex induced subgraphs on a
# small graph. Returns a per-class count table — used as ground truth.
# O(N^4); only used on tiny fixtures.
brute_size4 <- function(nodes, edges) {
  ids <- nodes$cell_id
  ekey <- paste(pmin(edges$source, edges$target),
                pmax(edges$source, edges$target), sep = "|")
  has_edge <- function(a, b) {
    paste(min(a, b), max(a, b), sep = "|") %in% ekey
  }
  classes <- character(0)
  N <- length(ids)
  for (i in seq_len(N - 3L)) {
    for (j in seq.int(i + 1L, N - 2L)) {
      for (k in seq.int(j + 1L, N - 1L)) {
        for (l in seq.int(k + 1L, N)) {
          v <- c(ids[i], ids[j], ids[k], ids[l])
          e <- c(
            has_edge(v[1], v[2]), has_edge(v[1], v[3]),
            has_edge(v[1], v[4]), has_edge(v[2], v[3]),
            has_edge(v[2], v[4]), has_edge(v[3], v[4])
          )
          ne <- sum(e)
          if (ne < 3L) next
          d <- c(
            e[1] + e[2] + e[3],   # v1
            e[1] + e[4] + e[5],   # v2
            e[2] + e[4] + e[6],   # v3
            e[3] + e[5] + e[6]    # v4
          )
          # Connected only if all degrees >= 1
          if (any(d == 0L)) next
          ds <- sort(d)
          cls <- if (ne == 3L && identical(ds, c(1L, 1L, 1L, 3L))) {
            "claw"
          } else if (ne == 3L && identical(ds, c(1L, 1L, 2L, 2L))) {
            "path"
          } else if (ne == 4L && identical(ds, c(2L, 2L, 2L, 2L))) {
            "cycle"
          } else if (ne == 4L && identical(ds, c(1L, 2L, 2L, 3L))) {
            "paw"
          } else if (ne == 5L && identical(ds, c(2L, 2L, 3L, 3L))) {
            "diamond"
          } else if (ne == 6L && identical(ds, c(3L, 3L, 3L, 3L))) {
            "K4"
          } else NA_character_
          if (!is.na(cls)) classes <- c(classes, cls)
        }
      }
    }
  }
  tab <- table(classes)
  data.table::data.table(
    canonical_iso_class = names(tab),
    count               = as.integer(tab)
  )[order(canonical_iso_class)]
}

# Build a SpatialGraph from a tiny edge list (helper to keep tests terse).
small_sg <- function(node_ids, edge_pairs, types = NULL) {
  if (is.null(types)) types <- "T"
  if (length(types) == 1L) types <- rep(types, length(node_ids))
  nodes <- data.frame(
    cell_id = node_ids,
    x = seq_along(node_ids), y = 0,
    cell_type = types,
    sample_id = "s1",
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    source = edge_pairs[, 1L], target = edge_pairs[, 2L],
    sample_id = "s1", stringsAsFactors = FALSE
  )
  SpatialGraph(nodes, edges)
}

test_that("each topology class is detected correctly on its canonical fixture", {
  # K4: 1 K4
  sg_k4 <- small_sg(letters[1:4], rbind(
    c("a","b"), c("a","c"), c("a","d"),
    c("b","c"), c("b","d"), c("c","d")
  ))
  cat_k4 <- find_motifs(sg_k4, size = 4L, colored = FALSE)$catalog
  expect_equal(cat_k4$canonical_iso_class, "K4")
  expect_equal(cat_k4$count, 1L)

  # P4 (path of 4): 1 path
  sg_p4 <- small_sg(letters[1:4], rbind(c("a","b"), c("b","c"), c("c","d")))
  cat_p4 <- find_motifs(sg_p4, size = 4L, colored = FALSE)$catalog
  expect_equal(cat_p4$canonical_iso_class, "path")
  expect_equal(cat_p4$count, 1L)

  # Claw: 1 claw
  sg_st <- small_sg(letters[1:4],
                     rbind(c("a","b"), c("a","c"), c("a","d")))
  cat_st <- find_motifs(sg_st, size = 4L, colored = FALSE)$catalog
  expect_equal(cat_st$canonical_iso_class, "claw")
  expect_equal(cat_st$count, 1L)

  # C4: 1 cycle
  sg_c4 <- small_sg(letters[1:4],
                     rbind(c("a","b"), c("b","c"), c("c","d"), c("a","d")))
  cat_c4 <- find_motifs(sg_c4, size = 4L, colored = FALSE)$catalog
  expect_equal(cat_c4$canonical_iso_class, "cycle")
  expect_equal(cat_c4$count, 1L)

  # Paw: triangle + pendant. 1 paw.
  sg_paw <- small_sg(letters[1:4],
                      rbind(c("a","b"), c("b","c"), c("a","c"), c("c","d")))
  cat_paw <- find_motifs(sg_paw, size = 4L, colored = FALSE)$catalog
  expect_equal(cat_paw$canonical_iso_class, "paw")
  expect_equal(cat_paw$count, 1L)

  # Diamond (K4 - one edge): 1 diamond.
  sg_dm <- small_sg(letters[1:4], rbind(
    c("a","b"), c("a","c"), c("a","d"),
    c("b","c"), c("b","d")
  ))
  cat_dm <- find_motifs(sg_dm, size = 4L, colored = FALSE)$catalog
  expect_equal(cat_dm$canonical_iso_class, "diamond")
  expect_equal(cat_dm$count, 1L)
})

test_that("known-count fixtures match expected counts per class", {
  # K_{2,3}: 3 cycles + 2 claws
  sg <- small_sg(
    letters[1:5],
    rbind(c("a","c"), c("a","d"), c("a","e"),
          c("b","c"), c("b","d"), c("b","e"))
  )
  m <- find_motifs(sg, size = 4L, colored = FALSE)
  cat_by_cls <- m$catalog[order(canonical_iso_class)]
  expect_equal(cat_by_cls$canonical_iso_class, c("claw", "cycle"))
  expect_equal(cat_by_cls$count, c(2L, 3L))

  # K5: 5 K4 instances (C(5,4)).
  pairs <- t(utils::combn(letters[1:5], 2L))
  sg_k5 <- small_sg(letters[1:5], pairs)
  cat_k5 <- find_motifs(sg_k5, size = 4L, colored = FALSE)$catalog
  expect_equal(cat_k5$canonical_iso_class, "K4")
  expect_equal(cat_k5$count, 5L)
})

test_that("brute-force vs backend on a small random graph (uncolored)", {
  fx <- fixture_random(N = 15L, seed = 13L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                             cell_id = fx$cell_id, method = "knn", k = 3L)
  m <- find_motifs(sg, size = 4L, colored = FALSE)
  bf <- brute_size4(sg$nodes, sg$edges)
  agg <- m$catalog[, .(count = sum(count)),
                   by = canonical_iso_class][order(canonical_iso_class)]
  expect_equal(agg$canonical_iso_class, bf$canonical_iso_class)
  expect_equal(agg$count, bf$count)
})

test_that("igraph::motifs(g, size = 4) total matches our enumerator total", {
  fx <- fixture_random(N = 25L, seed = 4L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                             cell_id = fx$cell_id, method = "knn", k = 4L)
  g <- as_igraph(sg)
  ig_counts <- igraph::motifs(g, size = 4L)
  ig_total <- sum(ig_counts, na.rm = TRUE)
  m <- find_motifs(sg, size = 4L, colored = FALSE)
  expect_equal(sum(m$catalog$count), ig_total)
})

test_that("colored size-4 motif total matches uncolored total", {
  fx <- fixture_random(N = 20L, seed = 9L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                             cell_id = fx$cell_id, method = "knn", k = 4L)
  m_un <- find_motifs(sg, size = 4L, colored = FALSE)
  m_col <- find_motifs(sg, size = 4L, colored = TRUE)
  expect_equal(sum(m_un$catalog$count), sum(m_col$catalog$count))
  # Colored partition agrees with uncolored when grouped by topology class.
  un_by <- m_un$catalog[, .(count = sum(count)),
                        by = canonical_iso_class][order(canonical_iso_class)]
  col_by <- m_col$catalog[, .(count = sum(count)),
                          by = canonical_iso_class][order(canonical_iso_class)]
  expect_equal(un_by, col_by)
})

test_that("anchored_on restricts size-4 enumeration", {
  fx <- fixture_random(N = 25L, seed = 6L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                             cell_id = fx$cell_id, method = "knn", k = 4L)
  anchors <- fx$cell_id[1:3]
  m <- find_motifs(sg, size = 4L, anchored_on = anchors, colored = FALSE)
  by_inst <- m$incidence[, .(touches = any(cell_id %in% anchors)),
                         by = instance_id]
  expect_true(all(by_inst$touches))
})

test_that("size-4 incidence has 4 rows per instance and orbits are consistent", {
  fx <- fixture_random(N = 20L, seed = 8L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                             cell_id = fx$cell_id, method = "knn", k = 4L)
  m <- find_motifs(sg, size = 4L, colored = TRUE)
  per_inst <- m$incidence[, .(n = .N), by = instance_id]
  expect_true(all(per_inst$n == 4L))
  per_inst_orb <- m$incidence[, .(sig = paste(sort(orbit_id),
                                              collapse = ",")),
                              by = .(motif_id, instance_id)]
  per_motif <- per_inst_orb[, .(uniq = data.table::uniqueN(sig)),
                            by = motif_id]
  expect_true(all(per_motif$uniq == 1L))
})

test_that("label-permutation null works for size 4 (uncolored topology total invariant)", {
  fx <- fixture_random(N = 25L, seed = 21L)
  sg <- build_spatial_graph(fx[, c("x","y")], fx$cell_type, fx$sample_id,
                             cell_id = fx$cell_id, method = "knn", k = 4L)
  m <- find_motifs(sg, size = 4L, colored = FALSE)
  res <- test_motif_enrichment(m, sg, n_perm = 50L, seed = 1L)
  # Label permutation cannot change topology counts when colored = FALSE.
  expect_equal(res$mean_null, res$observed)
  expect_equal(res$sd_null, rep(0, nrow(res)))
})
