#' Find motifs in a `SpatialGraph`
#'
#' Enumerate all motifs of a given `size` in `sg`, returning a `Motifs`
#' S3 object. v0.1 supports sizes 2 and 3 (size 4 lands in step 8). All
#' connected sub-graphs of that size are enumerated; the resulting catalog
#' distinguishes topology classes (e.g., closed triangle vs open wedge for
#' size 3) via the `canonical_iso_class` column.
#'
#' When `colored = TRUE` (default), motifs are further keyed by the sorted
#' tuple of cell-type colors per orbit, so e.g. a `T-T-B` triangle is a
#' different motif from a `T-B-B` triangle. The canonical labeling sorts
#' vertices by `(degree-in-motif, color)`, which is exact for sizes 2 and
#' 3 combined with topology class.
#'
#' @param sg a `SpatialGraph`.
#' @param size integer; one of 2, 3 (4 in v0.2).
#' @param backend `"igraph"` (only available in v0.1) or `"rust"` (planned).
#' @param colored logical; if `TRUE`, distinguish motifs by ordered color
#'   tuples in addition to topology.
#' @param anchored_on optional character vector of `cell_id`s. When set,
#'   only motifs containing at least one anchored cell are returned —
#'   intended to keep enumeration tractable on large graphs (size >= 4).
#' @param max_instances optional integer cap on enumerated instances per
#'   motif. If exceeded, a deterministic random sample of size `max_instances`
#'   is kept and a warning is emitted.
#' @return a `Motifs` object — a list with `catalog`, `incidence`, `meta`.
#' @examples
#' set.seed(1)
#' coords <- cbind(x = runif(20), y = runif(20))
#' sg <- build_spatial_graph(
#'   coords, cell_types = sample(c("T","B"), 20, TRUE),
#'   sample_id = "s1", method = "knn", k = 3
#' )
#' m <- find_motifs(sg, size = 2)
#' m$catalog
#' @export
find_motifs <- function(sg,
                        size = 3L,
                        backend = "igraph",
                        colored = TRUE,
                        anchored_on = NULL,
                        max_instances = NULL) {
  if (!inherits(sg, "SpatialGraph")) {
    stop("`sg` must be a SpatialGraph", call. = FALSE)
  }
  size <- as.integer(size)
  if (length(size) != 1L || is.na(size) || !size %in% c(2L, 3L, 4L)) {
    stop("`size` must be one of 2, 3, 4", call. = FALSE)
  }
  backend <- match.arg(backend, c("igraph", "rust"))
  if (identical(backend, "rust")) {
    if (!requireNamespace("smotifrs", quietly = TRUE)) {
      stop(
        "backend = 'rust' requires the smotifrs companion package. ",
        "Install pre-built binaries via R-universe (no Rust toolchain ",
        "needed): install.packages('smotifrs', ",
        "repos = c('https://rdries.r-universe.dev', getOption('repos')))",
        call. = FALSE
      )
    }
    return(smotifrs::find_motifs_rs(
      sg = sg, size = size, colored = colored,
      anchored_on = anchored_on, max_instances = max_instances
    ))
  }
  if (!is.null(anchored_on)) {
    anchored_on <- as.character(anchored_on)
    miss <- setdiff(anchored_on, sg$nodes$cell_id)
    if (length(miss)) {
      stop(length(miss),
           " cell_id(s) in anchored_on are not in sg$nodes$cell_id",
           call. = FALSE)
    }
  }
  if (!is.null(max_instances)) {
    max_instances <- as.integer(max_instances)
    if (is.na(max_instances) || max_instances < 1L) {
      stop("max_instances must be a positive integer", call. = FALSE)
    }
  }
  .find_motifs_igraph(
    sg = sg, size = size, colored = colored,
    anchored_on = anchored_on, max_instances = max_instances
  )
}

#' Construct a `Motifs` object
#'
#' Internal-style constructor used by `find_motifs()` to assemble its
#' return value. Kept exported because tests and downstream code may
#' build Motifs objects from custom enumerators.
#'
#' @param catalog a `data.table` with `motif_id, size, canonical_iso_class, color_tuple, count`.
#' @param incidence a `data.table` with `motif_id, instance_id, cell_id, orbit_id`.
#' @param instance_meta optional named list of per-topology-class wide
#'   instance tables, used by [test_motif_enrichment()] for fast
#'   permutation recounts. Expected names: `edge`, `triangle`, `wedge`.
#'   Each table holds the cell ids of an instance in canonical role
#'   order — e.g. `wedge` has `instance_id, end1, end2, center`. May be
#'   `NULL` for `Motifs` objects produced by hand.
#' @param meta a named list (`size`, `backend`, `colored`, `sg_hash`).
#' @return a `Motifs` S3 object.
#' @export
Motifs <- function(catalog, incidence, instance_meta = NULL,
                   meta = list()) {
  catalog <- data.table::as.data.table(catalog)
  incidence <- data.table::as.data.table(incidence)
  req_cat <- c("motif_id", "size", "canonical_iso_class",
               "color_tuple", "count")
  req_inc <- c("motif_id", "instance_id", "cell_id", "orbit_id")
  m_cat <- setdiff(req_cat, names(catalog))
  m_inc <- setdiff(req_inc, names(incidence))
  if (length(m_cat)) {
    stop("catalog missing columns: ", paste(m_cat, collapse = ", "),
         call. = FALSE)
  }
  if (length(m_inc)) {
    stop("incidence missing columns: ", paste(m_inc, collapse = ", "),
         call. = FALSE)
  }
  if (nrow(catalog)) data.table::setkey(catalog, motif_id)
  if (nrow(incidence)) data.table::setkey(incidence, motif_id, instance_id)
  if (is.null(instance_meta)) instance_meta <- list()
  structure(
    list(
      catalog       = catalog,
      incidence     = incidence,
      instance_meta = instance_meta,
      meta          = meta
    ),
    class = "Motifs"
  )
}

#' Bar chart of the most common motifs
#'
#' Returns a `ggplot2` object (requires the optional `ggplot2` Suggests).
#' Bars are sorted by count, descending. For colored motifs the labels
#' include the color tuple, which can get long — a default truncation
#' keeps display readable.
#'
#' @param x a `Motifs` object.
#' @param top integer; number of top motifs to plot.
#' @param max_label_chars integer; truncate `motif_id` labels longer
#'   than this with an ellipsis.
#' @param ... ignored.
#' @return a `ggplot` object.
#' @export
plot.Motifs <- function(x, top = 20L, max_label_chars = 60L, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot.Motifs(); install.packages('ggplot2')",
         call. = FALSE)
  }
  if (!nrow(x$catalog)) {
    stop("Motifs catalog is empty; nothing to plot", call. = FALSE)
  }
  dt <- data.table::copy(x$catalog)
  data.table::setorder(dt, -count)
  dt <- utils::head(dt, top)
  long <- nchar(dt$motif_id) > max_label_chars
  if (any(long)) {
    dt[long, motif_id := paste0(substr(motif_id, 1L, max_label_chars - 3L),
                                 "...")]
  }
  ggplot2::ggplot(
    dt,
    ggplot2::aes(x = stats::reorder(motif_id, count), y = count)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL, y = "instance count",
      title = sprintf("Top %d motifs (size %s)",
                      nrow(dt),
                      as.character(x$meta$size %||% NA))
    ) +
    ggplot2::theme_minimal()
}

#' @rdname Motifs
#' @param x a `Motifs` object.
#' @param ... ignored.
#' @export
print.Motifs <- function(x, ...) {
  cat("<Motifs>\n")
  cat(sprintf("  size: %d  |  colored: %s  |  backend: %s\n",
              x$meta$size %||% NA_integer_,
              if (isTRUE(x$meta$colored)) "yes" else "no",
              x$meta$backend %||% "?"))
  cat(sprintf("  motif types: %d  |  total instances: %d\n",
              nrow(x$catalog),
              if (nrow(x$catalog)) sum(x$catalog$count) else 0L))
  topn <- utils::head(x$catalog[order(-count)], 5L)
  if (nrow(topn)) {
    cat("  top motifs:\n")
    for (i in seq_len(nrow(topn))) {
      cat(sprintf("    %s  (count = %d)\n",
                  topn$motif_id[i], topn$count[i]))
    }
  }
  invisible(x)
}

# ---- internal: shared helpers for motif backends ---------------------------

# Compute a per-vertex orbit_id from a vector of (degree-in-motif, color)
# tuples. Vertices sharing a tuple share an orbit; orbit_id is the rank of
# the unique tuple in lexicographic order. Approximate canonical labeling,
# exact for sizes <= 4 combined with topology class — see motifs.R docs.
.orbit_from_degree_color <- function(degree, color) {
  key <- paste(degree, color, sep = "\x1f")
  uniq <- sort(unique(key))
  match(key, uniq)
}

# Hash of a SpatialGraph — used for cache invalidation in `meta$sg_hash`.
# Deliberately cheap: we hash node and edge identifiers, not coordinates.
.sg_hash <- function(sg) {
  .smotif_hash(list(
    sort(sg$nodes$cell_id),
    sg$edges[, paste0(source, "\x1f", target)]
  ))
}

# Cap the total number of motif instances (collectively, across motif_ids)
# at `max_instances`. Random subsample is seeded for reproducibility.
# `instance_meta` is filtered in lockstep so the per-class wide tables
# stay aligned with the surviving instances.
.cap_instances <- function(catalog, incidence, instance_meta, max_instances) {
  if (is.null(max_instances)) {
    return(list(catalog = catalog, incidence = incidence,
                instance_meta = instance_meta))
  }
  total <- if (nrow(catalog)) sum(catalog$count) else 0L
  if (total <= max_instances) {
    return(list(catalog = catalog, incidence = incidence,
                instance_meta = instance_meta))
  }
  warning(sprintf(
    "max_instances = %d is below total enumerated count (%d); subsampling",
    max_instances, total
  ), call. = FALSE)
  set.seed(1L)
  all_inst <- unique(incidence$instance_id)
  keep <- sort(sample(all_inst, max_instances, replace = FALSE))
  incidence <- incidence[instance_id %in% keep]
  catalog <- incidence[, .(count = data.table::uniqueN(instance_id)),
                       by = motif_id] |>
    merge(unique(catalog[, .(motif_id, size, canonical_iso_class,
                             color_tuple)]),
          by = "motif_id")
  data.table::setcolorder(catalog,
    c("motif_id", "size", "canonical_iso_class", "color_tuple", "count"))
  for (cls in names(instance_meta)) {
    instance_meta[[cls]] <- instance_meta[[cls]][instance_id %in% keep]
  }
  list(catalog = catalog, incidence = incidence,
       instance_meta = instance_meta)
}
