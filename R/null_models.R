# Null-model generators for motif enrichment. Three flavors are supported
# in v0.1:
#   * `label_permutation`  — shuffle `cell_type` across all nodes; the
#     graph topology is unchanged.
#   * `region_stratified`  — shuffle `cell_type` within each `region_id`;
#     preserves per-region cell-type composition.
#   * `geometric_jitter`   — perturb (x, y) by Gaussian noise with sd =
#     `null_args$sd`, rebuild the graph using the original construction
#     parameters, recount.
#
# `label_permutation` and `region_stratified` use the per-class wide
# instance tables on the `Motifs` object to recount in O(n_instances)
# without re-enumerating the graph. `geometric_jitter` is the expensive
# null — it re-runs `build_spatial_graph` and `find_motifs` per draw.

# Permute a cell_type vector globally.
.permute_labels_global <- function(cell_type) {
  cell_type[sample.int(length(cell_type))]
}

# Permute cell_type within each region. Cells with NA region_id are
# shuffled together (treated as their own implicit region). Preserves
# per-region cell-type composition exactly.
.permute_labels_stratified <- function(cell_type, region_id) {
  if (length(region_id) != length(cell_type)) {
    stop("cell_type and region_id must have equal length", call. = FALSE)
  }
  region_id <- as.character(region_id)
  region_id[is.na(region_id)] <- "<NA>"
  out <- cell_type
  for (rg in unique(region_id)) {
    idx <- which(region_id == rg)
    if (length(idx) > 1L) out[idx] <- cell_type[sample(idx)]
  }
  out
}

# Compute, for each topology class table, a vector of motif_ids under a
# given cell_id -> color lookup. Uses the per-class wide instance tables
# carried on the `Motifs` object (see motifs.R `instance_meta`).
.recount_labels <- function(instance_meta, color_lookup, colored) {
  parts <- list()
  if (!is.null(instance_meta$edge) && nrow(instance_meta$edge)) {
    e <- instance_meta$edge
    c1 <- color_lookup[e$v1]
    c2 <- color_lookup[e$v2]
    parts$edge <- if (colored) {
      paste0("size2_edge_", pmin(c1, c2), "-", pmax(c1, c2))
    } else {
      rep("size2_edge", length(c1))
    }
  }
  if (!is.null(instance_meta$triangle) && nrow(instance_meta$triangle)) {
    tri <- instance_meta$triangle
    c1 <- color_lookup[tri$v1]
    c2 <- color_lookup[tri$v2]
    c3 <- color_lookup[tri$v3]
    parts$triangle <- if (colored) {
      lvls <- sort(unique(c(c1, c2, c3)))
      k1 <- match(c1, lvls); k2 <- match(c2, lvls); k3 <- match(c3, lvls)
      k_lo <- pmin(k1, k2, k3)
      k_hi <- pmax(k1, k2, k3)
      k_md <- k1 + k2 + k3 - k_lo - k_hi
      paste0("size3_closed_",
             lvls[k_lo], "-", lvls[k_md], "-", lvls[k_hi])
    } else {
      rep("size3_closed", length(c1))
    }
  }
  if (!is.null(instance_meta$wedge) && nrow(instance_meta$wedge)) {
    w <- instance_meta$wedge
    c_e1 <- color_lookup[w$end1]
    c_e2 <- color_lookup[w$end2]
    c_ct <- color_lookup[w$center]
    parts$wedge <- if (colored) {
      paste0("size3_open_",
             pmin(c_e1, c_e2), "-", pmax(c_e1, c_e2), "_", c_ct)
    } else {
      rep("size3_open", length(c_e1))
    }
  }
  # Size-4 classes: each carries (w1..w4, d1..d4) so we can recompute
  # the (deg, color) sorted tuple under the new color assignment.
  for (cls in c("claw", "path", "cycle", "paw", "diamond", "K4")) {
    sub <- instance_meta[[cls]]
    if (is.null(sub) || !nrow(sub)) next
    if (colored) {
      cs <- .sort_colors_by_degree(
        sub$d1, sub$d2, sub$d3, sub$d4,
        color_lookup[sub$w1], color_lookup[sub$w2],
        color_lookup[sub$w3], color_lookup[sub$w4]
      )
      parts[[cls]] <- paste0(
        "size4_", cls, "_",
        cs$c1, "-", cs$c2, "-", cs$c3, "-", cs$c4
      )
    } else {
      parts[[cls]] <- rep(paste0("size4_", cls), nrow(sub))
    }
  }
  unlist(parts, use.names = FALSE)
}

# Recount per-motif counts under a permuted color assignment. Returns a
# named integer vector keyed on motif_id, restricted to the motif_ids
# present in `observed_motif_ids` (filling 0 for motifs that vanish).
.recount_one_perm <- function(instance_meta, color_lookup, colored,
                              observed_motif_ids) {
  ids <- .recount_labels(instance_meta, color_lookup, colored)
  tab <- table(ids)
  out <- rep(0L, length(observed_motif_ids))
  names(out) <- observed_motif_ids
  hit <- intersect(observed_motif_ids, names(tab))
  out[hit] <- as.integer(tab[hit])
  out
}

# Single label-shuffle null draw, returning a per-motif count vector.
.draw_label_perm <- function(motifs, sg, observed_motif_ids,
                             stratified, null_args) {
  ct <- sg$nodes$cell_type
  cell_ids <- sg$nodes$cell_id
  ct_perm <- if (stratified) {
    if (!"region_id" %in% names(sg$nodes)) {
      stop("region_stratified null requires a region_id column on sg$nodes",
           call. = FALSE)
    }
    .permute_labels_stratified(ct, sg$nodes$region_id)
  } else {
    .permute_labels_global(ct)
  }
  color_lookup <- ct_perm
  names(color_lookup) <- cell_ids
  .recount_one_perm(motifs$instance_meta, color_lookup,
                    isTRUE(motifs$meta$colored), observed_motif_ids)
}

# Single geometric-jitter null draw. Perturbs coords, rebuilds the graph
# with the same method/k/radius the original sg used, and re-runs
# `find_motifs` of the same size. Slow.
.draw_jitter <- function(motifs, sg, observed_motif_ids, null_args) {
  sd <- null_args$sd
  if (is.null(sd) || !is.finite(sd) || sd <= 0) {
    stop("geometric_jitter requires null_args$sd > 0 ", call. = FALSE)
  }
  build <- sg$meta$build_method
  if (is.null(build) || !build %in% c("knn", "delaunay", "radius")) {
    stop(
      "geometric_jitter needs sg to have been produced by ",
      "build_spatial_graph() so it can be rebuilt; ",
      "current build_method = '", build %||% "?",
      "'. Either rebuild the SpatialGraph or write a custom null.",
      call. = FALSE
    )
  }
  N <- nrow(sg$nodes)
  jx <- sg$nodes$x + stats::rnorm(N, sd = sd)
  jy <- sg$nodes$y + stats::rnorm(N, sd = sd)
  coords <- cbind(x = jx, y = jy)
  k_arg <- sg$meta$k
  if (is.null(k_arg) || is.na(k_arg)) k_arg <- 6L
  rad_arg <- sg$meta$radius
  if (is.null(rad_arg) || is.na(rad_arg)) rad_arg <- NULL
  md_arg <- sg$meta$max_distance
  if (is.null(md_arg) || is.na(md_arg)) md_arg <- NULL
  sg_new <- build_spatial_graph(
    coords      = coords,
    cell_types  = sg$nodes$cell_type,
    sample_id   = sg$nodes$sample_id,
    region_id   = sg$nodes$region_id,
    cell_id     = sg$nodes$cell_id,
    method      = build,
    k           = k_arg,
    radius      = rad_arg,
    max_distance = md_arg
  )
  m_new <- find_motifs(
    sg_new, size = motifs$meta$size, backend = "igraph",
    colored = isTRUE(motifs$meta$colored)
  )
  cnt <- m_new$catalog$count
  names(cnt) <- m_new$catalog$motif_id
  out <- rep(0L, length(observed_motif_ids))
  names(out) <- observed_motif_ids
  hit <- intersect(observed_motif_ids, names(cnt))
  out[hit] <- as.integer(cnt[hit])
  out
}
