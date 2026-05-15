# igraph-backed motif enumeration. Three topology classes are wired up
# in v0.1: edges (size 2), triangles (size 3 closed), and wedges/paths
# (size 3 open). Each backend returns (catalog, incidence) data.tables
# matching the schema expected by `Motifs()`. Size 4 lands in step 8.

# Top-level dispatcher invoked from `find_motifs()` once arguments have
# been validated.
.find_motifs_igraph <- function(sg, size, colored, anchored_on,
                                max_instances) {
  ct <- stats::setNames(sg$nodes$cell_type, sg$nodes$cell_id)

  parts <- if (size == 2L) {
    list(.find_motifs_size2(sg, ct, colored, anchored_on))
  } else if (size == 3L) {
    g <- as_igraph(sg)
    list(
      .find_motifs_size3_closed(sg, g, ct, colored, anchored_on),
      .find_motifs_size3_open(sg, ct, colored, anchored_on)
    )
  } else if (size == 4L) {
    g <- as_igraph(sg)
    list(.find_motifs_size4(sg, g, ct, colored, anchored_on))
  } else {
    stop("internal: unsupported size ", size, call. = FALSE)
  }

  catalog <- data.table::rbindlist(lapply(parts, `[[`, "catalog"),
                                   fill = TRUE)
  incidence <- data.table::rbindlist(lapply(parts, `[[`, "incidence"),
                                     fill = TRUE)
  instance_meta <- list()
  for (p in parts) {
    if (!is.null(p$instance_meta)) {
      instance_meta <- utils::modifyList(instance_meta, p$instance_meta)
    }
  }
  capped <- .cap_instances(catalog, incidence, instance_meta, max_instances)

  Motifs(
    catalog       = capped$catalog,
    incidence     = capped$incidence,
    instance_meta = capped$instance_meta,
    meta = list(
      size         = size,
      backend      = "igraph",
      colored      = colored,
      anchored     = !is.null(anchored_on),
      sg_hash      = .sg_hash(sg),
      generated_at = Sys.time()
    )
  )
}

# ---- shared helpers -------------------------------------------------------

# Assign orbit_id within each motif_id from per-vertex (degree, color)
# pairs. Approximate canonical labeling: vertices sharing a (degree,
# color) tuple share an orbit; orbit_id is the rank of the unique tuple
# in lexicographic order. Within a motif_id every instance has the same
# (degree, color) pattern by construction, so the orbit assignment is
# consistent.
.assign_orbit_ids <- function(inc, colored) {
  if (!colored) {
    # Without coloring, orbit_id collapses onto degree-in-motif alone.
    inc[, orbit_id := match(deg, sort(unique(deg))), by = motif_id]
  } else {
    inc[, orbit_id := .orbit_from_degree_color(deg, color),
        by = motif_id]
  }
  inc
}

# ---- size 2: edges --------------------------------------------------------

.find_motifs_size2 <- function(sg, ct, colored, anchored_on) {
  e <- data.table::copy(sg$edges)
  if (!nrow(e)) return(.empty_motif_part())

  if (!is.null(anchored_on)) {
    e <- e[source %in% anchored_on | target %in% anchored_on]
    if (!nrow(e)) return(.empty_motif_part())
  }

  e[, ca := ct[source]]
  e[, cb := ct[target]]

  if (colored) {
    e[, color_tuple := paste(pmin(ca, cb), pmax(ca, cb), sep = "-")]
    e[, motif_id := paste0("size2_edge_", color_tuple)]
  } else {
    e[, color_tuple := NA_character_]
    e[, motif_id := "size2_edge"]
  }

  data.table::setorder(e, motif_id, source, target)
  e[, instance_id := paste0("e_", .I)]

  catalog <- e[, .(count = .N), keyby = .(motif_id, color_tuple)]
  catalog[, `:=`(size = 2L, canonical_iso_class = "edge")]
  data.table::setcolorder(catalog,
    c("motif_id", "size", "canonical_iso_class", "color_tuple", "count"))

  # Two rows per edge instance.
  inc <- data.table::rbindlist(list(
    e[, .(motif_id, instance_id, cell_id = source, color = ca, deg = 1L)],
    e[, .(motif_id, instance_id, cell_id = target, color = cb, deg = 1L)]
  ))
  inc <- .assign_orbit_ids(inc, colored)

  inst_edges <- e[, .(instance_id, v1 = source, v2 = target)]

  list(
    catalog   = catalog,
    incidence = inc[, .(motif_id, instance_id, cell_id, orbit_id)],
    instance_meta = list(edge = inst_edges)
  )
}

# ---- size 3 closed: triangles --------------------------------------------

.find_motifs_size3_closed <- function(sg, g, ct, colored, anchored_on) {
  tri_idx <- igraph::triangles(g)
  if (length(tri_idx) == 0L) return(.empty_motif_part())

  V <- igraph::V(g)$name
  m <- matrix(V[as.integer(tri_idx)], ncol = 3L, byrow = TRUE,
              dimnames = list(NULL, c("v1", "v2", "v3")))

  if (!is.null(anchored_on)) {
    keep <- (m[, 1] %in% anchored_on |
             m[, 2] %in% anchored_on |
             m[, 3] %in% anchored_on)
    m <- m[keep, , drop = FALSE]
    if (!nrow(m)) return(.empty_motif_part())
  }

  c1 <- ct[m[, 1]]; c2 <- ct[m[, 2]]; c3 <- ct[m[, 3]]
  # Sort the three colors per row in a vectorized way: encode as integer
  # codes against a shared factor level set, take pmin/pmid/pmax, decode
  # back to strings. Avoids a per-row apply.
  lvls <- sort(unique(c(c1, c2, c3)))
  k1 <- match(c1, lvls); k2 <- match(c2, lvls); k3 <- match(c3, lvls)
  k_lo <- pmin(k1, k2, k3)
  k_hi <- pmax(k1, k2, k3)
  k_md <- k1 + k2 + k3 - k_lo - k_hi
  cols_sorted <- cbind(lvls[k_lo], lvls[k_md], lvls[k_hi])

  if (colored) {
    color_tuple <- paste(cols_sorted[, 1], cols_sorted[, 2],
                         cols_sorted[, 3], sep = "-")
    motif_id <- paste0("size3_closed_", color_tuple)
  } else {
    color_tuple <- rep(NA_character_, nrow(m))
    motif_id <- rep("size3_closed", nrow(m))
  }

  inst_dt <- data.table::data.table(
    motif_id    = motif_id,
    color_tuple = color_tuple,
    v1 = m[, 1], v2 = m[, 2], v3 = m[, 3],
    c1 = c1, c2 = c2, c3 = c3
  )
  data.table::setorder(inst_dt, motif_id, v1, v2, v3)
  inst_dt[, instance_id := paste0("tri_", .I)]

  catalog <- inst_dt[, .(count = .N), keyby = .(motif_id, color_tuple)]
  catalog[, `:=`(size = 3L, canonical_iso_class = "triangle")]
  data.table::setcolorder(catalog,
    c("motif_id", "size", "canonical_iso_class", "color_tuple", "count"))

  # Triangle: all three vertices have degree 2 in the motif.
  inc <- data.table::rbindlist(list(
    inst_dt[, .(motif_id, instance_id, cell_id = v1, color = c1, deg = 2L)],
    inst_dt[, .(motif_id, instance_id, cell_id = v2, color = c2, deg = 2L)],
    inst_dt[, .(motif_id, instance_id, cell_id = v3, color = c3, deg = 2L)]
  ))
  inc <- .assign_orbit_ids(inc, colored)

  inst_tri <- inst_dt[, .(instance_id, v1, v2, v3)]

  list(
    catalog   = catalog,
    incidence = inc[, .(motif_id, instance_id, cell_id, orbit_id)],
    instance_meta = list(triangle = inst_tri)
  )
}

# ---- size 3 open: wedges (paths a–b–c with no a–c edge) -------------------

.find_motifs_size3_open <- function(sg, ct, colored, anchored_on) {
  if (!nrow(sg$edges)) return(.empty_motif_part())

  # Directed adjacency view: each undirected edge contributes two rows
  # (b -> x). Self-join on b yields all (a, b, c) triples sharing center b.
  adj <- data.table::rbindlist(list(
    sg$edges[, .(b = source, x = target)],
    sg$edges[, .(b = target, x = source)]
  ))
  data.table::setkey(adj, b)
  triples <- merge(adj, adj, by = "b", allow.cartesian = TRUE,
                   suffixes = c("_a", "_c"))
  # Enforce x_a < x_c so each unordered {a, c} pair appears once and
  # x_a == x_c (the same edge taken twice) is ruled out.
  triples <- triples[x_a < x_c]
  if (!nrow(triples)) return(.empty_motif_part())

  # Drop closed triangles: keep only triples where (a, c) is NOT an edge.
  edge_keys <- sg$edges[, paste(source, target, sep = "\x1f")]
  triples[, ac_key := paste(pmin(x_a, x_c), pmax(x_a, x_c),
                            sep = "\x1f")]
  triples <- triples[!ac_key %in% edge_keys]
  if (!nrow(triples)) return(.empty_motif_part())

  if (!is.null(anchored_on)) {
    keep <- (triples$x_a %in% anchored_on |
             triples$b %in% anchored_on |
             triples$x_c %in% anchored_on)
    triples <- triples[keep]
    if (!nrow(triples)) return(.empty_motif_part())
  }

  triples[, ca := ct[x_a]]
  triples[, cb := ct[b]]
  triples[, cc := ct[x_c]]
  # Canonical labeling: ends share degree 1, center has degree 2. Sort the
  # two end colors so the printed tuple is deterministic.
  triples[, end_lo := pmin(ca, cc)]
  triples[, end_hi := pmax(ca, cc)]

  if (colored) {
    triples[, color_tuple := paste0(end_lo, "-", end_hi, "_", cb)]
    triples[, motif_id := paste0("size3_open_", color_tuple)]
  } else {
    triples[, color_tuple := NA_character_]
    triples[, motif_id := "size3_open"]
  }

  data.table::setorder(triples, motif_id, b, x_a, x_c)
  triples[, instance_id := paste0("wedge_", .I)]

  catalog <- triples[, .(count = .N), keyby = .(motif_id, color_tuple)]
  catalog[, `:=`(size = 3L, canonical_iso_class = "wedge")]
  data.table::setcolorder(catalog,
    c("motif_id", "size", "canonical_iso_class", "color_tuple", "count"))

  # Three vertices per wedge: end-a (deg 1), center (deg 2), end-c (deg 1).
  inc <- data.table::rbindlist(list(
    triples[, .(motif_id, instance_id, cell_id = x_a, color = ca, deg = 1L)],
    triples[, .(motif_id, instance_id, cell_id = b,   color = cb, deg = 2L)],
    triples[, .(motif_id, instance_id, cell_id = x_c, color = cc, deg = 1L)]
  ))
  inc <- .assign_orbit_ids(inc, colored)

  inst_wedge <- triples[, .(instance_id, end1 = x_a, end2 = x_c, center = b)]

  list(
    catalog   = catalog,
    incidence = inc[, .(motif_id, instance_id, cell_id, orbit_id)],
    instance_meta = list(wedge = inst_wedge)
  )
}

# ---- size 4: connected 4-vertex subgraphs --------------------------------
#
# Strategy: for each connected 3-subgraph (triangle or wedge), extend by a
# fourth vertex d in the union of neighborhoods, dedupe via per-row sort
# of the 4 cell ids, then classify the induced subgraph by (n_edges,
# sorted degree sequence) into one of six topology classes:
#   * claw    (star K_{1,3})        E = 3, deg = (1,1,1,3)
#   * path    (P4)                  E = 3, deg = (1,1,2,2)
#   * cycle   (C4)                  E = 4, deg = (2,2,2,2)
#   * paw     (triangle + pendant)  E = 4, deg = (1,2,2,3)
#   * diamond (K4 minus one edge)   E = 5, deg = (2,2,3,3)
#   * K4      (complete)            E = 6, deg = (3,3,3,3)
#
# `igraph::motifs(g, size = 4)` is used in tests to cross-check the
# total count per topology class — the design treats it as ground truth
# for counts and our enumerator as the instance-level provider.

.find_motifs_size4 <- function(sg, g, ct, colored, anchored_on) {
  if (!nrow(sg$edges)) return(.empty_motif_part())

  # 1. Enumerate connected 3-subgraphs as canonical (v1 < v2 < v3).
  three_dt <- .enumerate_connected_3subgraphs(sg, g)
  if (!nrow(three_dt)) return(.empty_motif_part())

  # 2. Extend by every neighbor of each of the three vertices that
  #    isn't already in the 3-set.
  adj <- data.table::rbindlist(list(
    sg$edges[, .(b = source, x = target)],
    sg$edges[, .(b = target, x = source)]
  ))
  data.table::setkey(adj, b)
  exp_dt <- data.table::rbindlist(list(
    merge(three_dt, adj, by.x = "v1", by.y = "b",
          allow.cartesian = TRUE)[, .(v1, v2, v3, d = x)],
    merge(three_dt, adj, by.x = "v2", by.y = "b",
          allow.cartesian = TRUE)[, .(v1, v2, v3, d = x)],
    merge(three_dt, adj, by.x = "v3", by.y = "b",
          allow.cartesian = TRUE)[, .(v1, v2, v3, d = x)]
  ))
  exp_dt <- exp_dt[d != v1 & d != v2 & d != v3]
  if (!nrow(exp_dt)) return(.empty_motif_part())

  # 3. Sort the four cell ids per row to canonical (w1 < w2 < w3 < w4)
  #    using a vectorized 4-element sorting network.
  s4 <- .sort4_keys(exp_dt$v1, exp_dt$v2, exp_dt$v3, exp_dt$d)
  exp_dt[, `:=`(w1 = s4$v1, w2 = s4$v2, w3 = s4$v3, w4 = s4$v4)]
  uniq <- unique(exp_dt[, .(w1, w2, w3, w4)])

  # 4. Classify each unique 4-set by (n_edges, sorted degree sequence).
  edge_set <- sg$edges[, paste(source, target, sep = "\x1f")]
  hash <- function(a, b) paste(a, b, sep = "\x1f")
  uniq[, e12 := hash(w1, w2) %in% edge_set]
  uniq[, e13 := hash(w1, w3) %in% edge_set]
  uniq[, e14 := hash(w1, w4) %in% edge_set]
  uniq[, e23 := hash(w2, w3) %in% edge_set]
  uniq[, e24 := hash(w2, w4) %in% edge_set]
  uniq[, e34 := hash(w3, w4) %in% edge_set]
  uniq[, n_edges := as.integer(e12 + e13 + e14 + e23 + e24 + e34)]
  uniq[, d1 := as.integer(e12 + e13 + e14)]
  uniq[, d2 := as.integer(e12 + e23 + e24)]
  uniq[, d3 := as.integer(e13 + e23 + e34)]
  uniq[, d4 := as.integer(e14 + e24 + e34)]
  ds <- .sort4_keys(uniq$d1, uniq$d2, uniq$d3, uniq$d4)
  uniq[, `:=`(ds1 = ds$v1, ds2 = ds$v2, ds3 = ds$v3, ds4 = ds$v4)]

  uniq[, canonical_iso_class := data.table::fcase(
    n_edges == 3L & ds1 == 1L & ds2 == 1L & ds3 == 1L & ds4 == 3L, "claw",
    n_edges == 3L & ds1 == 1L & ds2 == 1L & ds3 == 2L & ds4 == 2L, "path",
    n_edges == 4L & ds1 == 2L & ds2 == 2L & ds3 == 2L & ds4 == 2L, "cycle",
    n_edges == 4L & ds1 == 1L & ds2 == 2L & ds3 == 2L & ds4 == 3L, "paw",
    n_edges == 5L & ds1 == 2L & ds2 == 2L & ds3 == 3L & ds4 == 3L, "diamond",
    n_edges == 6L & ds1 == 3L & ds2 == 3L & ds3 == 3L & ds4 == 3L, "K4"
  )]
  uniq <- uniq[!is.na(canonical_iso_class)]
  if (!nrow(uniq)) return(.empty_motif_part())

  if (!is.null(anchored_on)) {
    keep <- (uniq$w1 %in% anchored_on |
             uniq$w2 %in% anchored_on |
             uniq$w3 %in% anchored_on |
             uniq$w4 %in% anchored_on)
    uniq <- uniq[keep]
    if (!nrow(uniq)) return(.empty_motif_part())
  }

  # 5. Color tuples via approximate canonical labeling: sort the four
  #    vertices by (degree-in-motif, color), then concatenate the colors.
  uniq[, c1 := ct[w1]]
  uniq[, c2 := ct[w2]]
  uniq[, c3 := ct[w3]]
  uniq[, c4 := ct[w4]]

  if (colored) {
    sorted_colors <- .sort_colors_by_degree(
      uniq$d1, uniq$d2, uniq$d3, uniq$d4,
      uniq$c1, uniq$c2, uniq$c3, uniq$c4
    )
    uniq[, color_tuple := paste(sorted_colors$c1, sorted_colors$c2,
                                sorted_colors$c3, sorted_colors$c4,
                                sep = "-")]
    uniq[, motif_id := paste0("size4_", canonical_iso_class, "_",
                              color_tuple)]
  } else {
    uniq[, color_tuple := NA_character_]
    uniq[, motif_id := paste0("size4_", canonical_iso_class)]
  }

  data.table::setorder(uniq, motif_id, w1, w2, w3, w4)
  uniq[, instance_id := paste0("q_", .I)]

  catalog <- uniq[, .(count = .N,
                      canonical_iso_class = canonical_iso_class[1L]),
                  keyby = .(motif_id, color_tuple)]
  catalog[, size := 4L]
  data.table::setcolorder(catalog,
    c("motif_id", "size", "canonical_iso_class", "color_tuple", "count"))

  inc <- data.table::rbindlist(list(
    uniq[, .(motif_id, instance_id, cell_id = w1, color = c1, deg = d1)],
    uniq[, .(motif_id, instance_id, cell_id = w2, color = c2, deg = d2)],
    uniq[, .(motif_id, instance_id, cell_id = w3, color = c3, deg = d3)],
    uniq[, .(motif_id, instance_id, cell_id = w4, color = c4, deg = d4)]
  ))
  inc <- .assign_orbit_ids(inc, colored)

  # Per-class wide instance tables. Carry per-position (deg) so the
  # null-model fast path can recompute color tuples without re-sorting
  # by topology class.
  inst_meta <- list()
  for (cls in c("claw", "path", "cycle", "paw", "diamond", "K4")) {
    sub <- uniq[canonical_iso_class == cls]
    if (nrow(sub)) {
      inst_meta[[cls]] <- sub[, .(instance_id,
                                  w1, w2, w3, w4,
                                  d1, d2, d3, d4)]
    }
  }

  list(
    catalog       = catalog,
    incidence     = inc[, .(motif_id, instance_id, cell_id, orbit_id)],
    instance_meta = inst_meta
  )
}

# Enumerate every connected 3-vertex induced subgraph (triangle or
# wedge) as a sorted (v1 < v2 < v3) tuple. Used as the seed set for
# size-4 expansion.
.enumerate_connected_3subgraphs <- function(sg, g) {
  tri_idx <- igraph::triangles(g)
  if (length(tri_idx)) {
    V <- igraph::V(g)$name
    tri_m <- matrix(V[as.integer(tri_idx)], ncol = 3L, byrow = TRUE)
    s3 <- .sort3_keys(tri_m[, 1L], tri_m[, 2L], tri_m[, 3L])
    tri_dt <- data.table::data.table(v1 = s3$v1, v2 = s3$v2, v3 = s3$v3)
  } else {
    tri_dt <- data.table::data.table(
      v1 = character(), v2 = character(), v3 = character()
    )
  }

  adj <- data.table::rbindlist(list(
    sg$edges[, .(b = source, x = target)],
    sg$edges[, .(b = target, x = source)]
  ))
  data.table::setkey(adj, b)
  triples <- merge(adj, adj, by = "b", allow.cartesian = TRUE,
                   suffixes = c("_a", "_c"))
  triples <- triples[x_a < x_c]
  edge_keys <- sg$edges[, paste(source, target, sep = "\x1f")]
  triples[, ac_key := paste(pmin(x_a, x_c), pmax(x_a, x_c),
                            sep = "\x1f")]
  triples <- triples[!ac_key %in% edge_keys]
  if (!nrow(triples)) {
    return(unique(tri_dt))
  }
  s3 <- .sort3_keys(triples$x_a, triples$b, triples$x_c)
  wedge_dt <- data.table::data.table(v1 = s3$v1, v2 = s3$v2, v3 = s3$v3)

  unique(data.table::rbindlist(list(tri_dt, wedge_dt)))
}

# Vectorized 3-element sort via a tiny sorting network.
.sort3_keys <- function(a, b, c) {
  ab_lo <- pmin(a, b); ab_hi <- pmax(a, b)
  v1 <- pmin(ab_lo, c); v3a <- pmax(ab_lo, c)
  v2 <- pmin(ab_hi, v3a); v3 <- pmax(ab_hi, v3a)
  list(v1 = v1, v2 = v2, v3 = v3)
}

# Vectorized 4-element sort via the 5-comparison Bose-Nelson network:
# (1,2),(3,4),(1,3),(2,4),(2,3). Works for character or numeric input.
.sort4_keys <- function(a, b, c, d) {
  ab_lo <- pmin(a, b); ab_hi <- pmax(a, b)
  cd_lo <- pmin(c, d); cd_hi <- pmax(c, d)
  v1 <- pmin(ab_lo, cd_lo); v3a <- pmax(ab_lo, cd_lo)
  v2a <- pmin(ab_hi, cd_hi); v4 <- pmax(ab_hi, cd_hi)
  v2 <- pmin(v2a, v3a); v3 <- pmax(v2a, v3a)
  list(v1 = v1, v2 = v2, v3 = v3, v4 = v4)
}

# Sort four (degree, color) pairs per row by the lexicographic
# (degree, color) order, returning the colors in sorted order. Degrees
# are zero-padded so the string comparison agrees with numeric ordering
# even when colors include digits.
.sort_colors_by_degree <- function(d1, d2, d3, d4, c1, c2, c3, c4) {
  k1 <- paste(sprintf("%02d", d1), c1, sep = "\x1f")
  k2 <- paste(sprintf("%02d", d2), c2, sep = "\x1f")
  k3 <- paste(sprintf("%02d", d3), c3, sep = "\x1f")
  k4 <- paste(sprintf("%02d", d4), c4, sep = "\x1f")
  s4 <- .sort4_keys(k1, k2, k3, k4)
  list(
    c1 = sub("^[^\x1f]+\x1f", "", s4$v1),
    c2 = sub("^[^\x1f]+\x1f", "", s4$v2),
    c3 = sub("^[^\x1f]+\x1f", "", s4$v3),
    c4 = sub("^[^\x1f]+\x1f", "", s4$v4)
  )
}

# ---- internal: empty result with the right schema ------------------------

.empty_motif_part <- function() {
  list(
    catalog = data.table::data.table(
      motif_id            = character(),
      size                = integer(),
      canonical_iso_class = character(),
      color_tuple         = character(),
      count               = integer()
    ),
    incidence = data.table::data.table(
      motif_id    = character(),
      instance_id = character(),
      cell_id     = character(),
      orbit_id    = integer()
    )
  )
}
