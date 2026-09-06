# motif_enrichment(): the enumerate-once / recount-many entry point.
#
# Prefers the Rust backend, which enumerates and tests in one call and never
# ships per-instance data back to R. Without it, a pure-R path covers sizes 2
# and 3: enumeration through igraph, canonical forms hard-coded for the three
# topologies that exist at those sizes, and the permutation recount done with
# tabulate() over integer class codes rather than by rebuilding a motif_id
# string per instance per draw. Size 4 in pure R is not offered -- it is one to
# two orders of magnitude more instances, and the honest answer is to install
# the backend rather than to run something that will not finish.

#' Motif enrichment on a spatial graph
#'
#' Enumerates every connected induced subgraph of the requested size and tests
#' each colored motif class against a permutation null. Enumeration happens
#' once; each permutation is a relabel-and-count over the fixed instance set.
#'
#' @param sg a [SpatialGraph()].
#' @param size motif size: 2, 3 or 4.
#' @param n_perm number of null draws.
#' @param seed integer seed.
#' @param null which null to test against. `"label"` shuffles labels over all
#'   nodes; `"stratified"` shuffles within `strata_column`; `"conditional"`
#'   holds the observed pairwise edge-type composition fixed, so a size-3 or
#'   size-4 motif that is still enriched is enriched *beyond* what its
#'   constituent pairs explain. See [smotifrs::motif_enrichment_rs()] for the
#'   sampler and its diagnostics.
#' @param strata_column name of a column on `sg$nodes` to stratify by.
#'   Defaults to `region_id` when that exists and `null = "stratified"`.
#' @param anchored_on optional character vector of `cell_id`s; enumeration is
#'   restricted to subgraphs touching at least one of them.
#' @param backend `"auto"`, `"rust"` or `"r"`.
#' @param ... passed to the backend (e.g. `cond_temp`).
#' @returns a `data.table`, one row per motif class: `motif_id`, `topology`,
#'   `size`, `color_tuple`, `observed`, `expected`, `sd_null`, `z`, `fold`,
#'   `p_enrich`, `p_deplete`, `p_adj`.
#' @examples
#' set.seed(1)
#' coords <- cbind(x = runif(60), y = runif(60))
#' sg <- build_spatial_graph(coords,
#'     cell_types = sample(c("T", "B"), 60, TRUE),
#'     sample_id = "s1", method = "knn", k = 4
#' )
#' motif_enrichment(sg, size = 3L, n_perm = 99L)
#' @export
motif_enrichment <- function(sg,
                             size = 3L,
                             n_perm = 1000L,
                             seed = 1L,
                             null = c("label", "stratified", "conditional"),
                             strata_column = NULL,
                             anchored_on = NULL,
                             backend = c("auto", "rust", "r"),
                             ...) {
  if (!inherits(sg, "SpatialGraph")) {
    stop("`sg` must be a SpatialGraph", call. = FALSE)
  }
  null <- match.arg(null)
  backend <- match.arg(backend)
  size <- as.integer(size)
  if (length(size) != 1L || is.na(size) || !size %in% 2:4) {
    stop("size must be 2, 3 or 4", call. = FALSE)
  }

  have_rs <- requireNamespace("smotifrs", quietly = TRUE)
  use_rust <- switch(backend,
    auto = have_rs,
    rust = TRUE,
    r = FALSE
  )
  if (backend == "rust" && !have_rs) {
    stop("backend = \"rust\" needs the smotifrs package", call. = FALSE)
  }
  if (!use_rust && size == 4L) {
    stop(
      "size 4 needs the smotifrs backend: it is one to two orders of ",
      "magnitude more instances than size 3, and the pure-R path would not ",
      "finish on a real dataset. Install smotifrs, or use size 2 or 3.",
      call. = FALSE
    )
  }
  if (!use_rust && null == "conditional") {
    stop(
      "the conditional null is implemented in the smotifrs backend only",
      call. = FALSE
    )
  }

  nodes <- sg$nodes
  ids <- as.character(nodes$cell_id)
  from <- match(as.character(sg$edges$source), ids)
  to <- match(as.character(sg$edges$target), ids)
  ct <- factor(as.character(nodes$cell_type))

  strata <- NULL
  if (null == "stratified") {
    scol <- strata_column %||% "region_id"
    if (!scol %in% names(nodes)) {
      stop(sprintf(
        "null = \"stratified\" needs column '%s' on sg$nodes", scol
      ), call. = FALSE)
    }
    strata <- factor(as.character(nodes[[scol]]))
  } else if (!is.null(strata_column)) {
    if (!strata_column %in% names(nodes)) {
      stop(sprintf("column '%s' not found on sg$nodes", strata_column),
        call. = FALSE
      )
    }
    strata <- factor(as.character(nodes[[strata_column]]))
    null <- "stratified"
  }

  anchors <- NULL
  if (!is.null(anchored_on)) {
    anchors <- match(as.character(anchored_on), ids)
    if (anyNA(anchors)) {
      stop("anchored_on contains cell_ids not present in sg$nodes",
        call. = FALSE
      )
    }
  }

  if (use_rust) {
    return(smotifrs::motif_enrichment_rs(
      from = from, to = to, cell_type = ct, n_nodes = length(ids),
      size = size, n_perm = n_perm, seed = seed, null = null,
      strata = strata, anchored_on = anchors, ...
    ))
  }
  .motif_enrichment_r(
    from = from, to = to, ct = ct, n_nodes = length(ids), size = size,
    n_perm = as.integer(n_perm), seed = as.integer(seed), strata = strata,
    anchors = anchors
  )
}


# --- pure-R path, sizes 2 and 3 ------------------------------------------
#
# Class codes are integers throughout. At size 3 there are exactly two
# topologies and their canonical forms are simple enough to write down:
#   open   (wedge)   -> the centre is distinguished; ends are interchangeable
#   closed (triangle)-> all three positions interchangeable, so sort
# which is what .cls3_* below encode. This mirrors what the Rust path derives
# from the automorphism group; the cross-backend test pins them together.

.enum_size2 <- function(ig) {
  e <- igraph::as_edgelist(ig, names = FALSE)
  list(a = e[, 1L], b = e[, 2L])
}

.enum_size3 <- function(ig) {
  tri <- matrix(as.integer(igraph::triangles(ig)), nrow = 3L)
  n <- igraph::vcount(ig)
  el <- igraph::as_edgelist(ig, names = FALSE)
  # Packed undirected edge keys for a vectorized membership test.
  # igraph::are_adjacent() takes single vertices, not vectors -- calling it on
  # a vector silently uses only the first element, which is what made an
  # earlier version of this undercount wedges by 3x.
  ekey <- function(a, b) pmin.int(a, b) * (n + 1) + pmax.int(a, b)
  edge_set <- ekey(el[, 1L], el[, 2L])
  adj <- igraph::as_adj_list(ig, mode = "all")
  cen <- vector("list", n)
  e1 <- vector("list", n)
  e2 <- vector("list", n)
  for (v in seq_len(n)) {
    nb <- as.integer(adj[[v]])
    if (length(nb) < 2L) next
    cb <- utils::combn(nb, 2L)
    keep <- !(ekey(cb[1L, ], cb[2L, ]) %in% edge_set)
    if (!any(keep)) next
    cen[[v]] <- rep.int(v, sum(keep))
    e1[[v]] <- cb[1L, keep]
    e2[[v]] <- cb[2L, keep]
  }
  list(
    tri = tri,
    w_center = unlist(cen, use.names = FALSE),
    w_e1 = unlist(e1, use.names = FALSE),
    w_e2 = unlist(e2, use.names = FALSE)
  )
}

.motif_enrichment_r <- function(from, to, ct, n_nodes, size, n_perm, seed,
                                strata, anchors) {
  K <- nlevels(ct)
  codes <- as.integer(ct)
  ig <- igraph::graph_from_edgelist(
    cbind(from, to),
    directed = FALSE
  )
  ig <- igraph::simplify(ig)
  if (igraph::vcount(ig) < n_nodes) {
    ig <- igraph::add_vertices(ig, n_nodes - igraph::vcount(ig))
  }

  if (size == 2L) {
    inst <- .enum_size2(ig)
    if (!is.null(anchors)) {
      keep <- inst$a %in% anchors | inst$b %in% anchors
      inst <- list(a = inst$a[keep], b = inst$b[keep])
    }
    nbins <- K * K
    tally <- function(cc) {
      lo <- pmin.int(cc[inst$a], cc[inst$b])
      hi <- pmax.int(cc[inst$a], cc[inst$b])
      tabulate((lo - 1L) * K + hi, nbins = nbins)
    }
    n_inst <- length(inst$a)
  } else {
    inst <- .enum_size3(ig)
    if (!is.null(anchors)) {
      if (ncol(inst$tri)) {
        kt <- apply(inst$tri, 2L, function(z) any(z %in% anchors))
        inst$tri <- inst$tri[, kt, drop = FALSE]
      }
      kw <- inst$w_center %in% anchors | inst$w_e1 %in% anchors |
        inst$w_e2 %in% anchors
      inst$w_center <- inst$w_center[kw]
      inst$w_e1 <- inst$w_e1[kw]
      inst$w_e2 <- inst$w_e2[kw]
    }
    # closed block first, then open; both use a K^3 address space
    blk <- K * K * K
    nbins <- 2L * blk
    tri <- inst$tri
    tally <- function(cc) {
      out <- integer(nbins)
      if (ncol(tri)) {
        c1 <- cc[tri[1L, ]]
        c2 <- cc[tri[2L, ]]
        c3 <- cc[tri[3L, ]]
        lo <- pmin.int(c1, c2, c3)
        hi <- pmax.int(c1, c2, c3)
        md <- c1 + c2 + c3 - lo - hi
        out <- out + tabulate(
          (lo - 1L) * K * K + (md - 1L) * K + hi,
          nbins = nbins
        )
      }
      if (length(inst$w_center)) {
        cc_c <- cc[inst$w_center]
        a <- pmin.int(cc[inst$w_e1], cc[inst$w_e2])
        b <- pmax.int(cc[inst$w_e1], cc[inst$w_e2])
        out <- out + tabulate(
          blk + (cc_c - 1L) * K * K + (a - 1L) * K + b,
          nbins = nbins
        )
      }
      out
    }
    n_inst <- ncol(tri) + length(inst$w_center)
  }

  obs <- tally(codes)
  sim <- matrix(0L, nrow = length(obs), ncol = n_perm)
  set.seed(seed)
  grp <- if (is.null(strata)) NULL else split(seq_len(n_nodes), strata)
  for (i in seq_len(n_perm)) {
    cc <- if (is.null(grp)) {
      codes[sample.int(n_nodes)]
    } else {
      z <- codes
      for (g in grp) if (length(g) > 1L) z[g] <- codes[sample(g)]
      z
    }
    sim[, i] <- tally(cc)
  }

  keep <- which(obs > 0L | rowSums(sim) > 0L)
  obs <- obs[keep]
  sim <- sim[keep, , drop = FALSE]
  lv <- levels(ct)

  if (size == 2L) {
    lo <- ((keep - 1L) %/% K) + 1L
    hi <- keep - (lo - 1L) * K
    topo <- rep("edge", length(keep))
    cols <- cbind(lv[lo], lv[hi])
  } else {
    is_open <- keep > blk
    idx <- ifelse(is_open, keep - blk, keep) - 1L
    i1 <- idx %/% (K * K) + 1L
    i2 <- (idx %% (K * K)) %/% K + 1L
    i3 <- idx %% K + 1L
    topo <- ifelse(is_open, "open", "closed")
    cols <- cbind(lv[i1], lv[i2], lv[i3])
  }

  mean_null <- rowMeans(sim)
  sd_null <- apply(sim, 1L, stats::sd)
  eps <- .Machine$double.eps
  out <- data.table::data.table(
    motif_id = paste0(
      "size", size, "_", topo, "_",
      apply(cols, 1L, paste, collapse = "-")
    ),
    topology = topo,
    size = size,
    color_tuple = split(cols, row(cols)),
    observed = as.numeric(obs),
    expected = mean_null,
    sd_null = sd_null,
    p_enrich = (1 + rowSums(sim >= obs)) / (1 + n_perm),
    p_deplete = (1 + rowSums(sim <= obs)) / (1 + n_perm)
  )
  names(out$color_tuple) <- NULL
  out[, "z" := (observed - expected) / pmax(sd_null, sqrt(eps))]
  out[, "fold" := observed / pmax(expected, eps)]
  out[, "p_adj" := stats::p.adjust(
    pmin(1, 2 * pmin(p_enrich, p_deplete)),
    method = "BH"
  )]
  data.table::setcolorder(out, c(
    "motif_id", "topology", "size", "color_tuple", "observed", "expected",
    "sd_null", "z", "fold", "p_enrich", "p_deplete", "p_adj"
  ))
  data.table::setorder(out, p_adj, -z)
  # double, not integer: instance counts pass 2^31 at size 4 on a
  # million-cell graph, and the Rust path returns a double too
  attr(out, "n_instances") <- as.numeric(n_inst)
  attr(out, "n_perm") <- n_perm
  attr(out, "null") <- if (is.null(strata)) "label" else "stratified"
  attr(out, "backend") <- "r"
  out[]
}
