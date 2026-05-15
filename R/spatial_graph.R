#' Spatial graph S3 class
#'
#' A `SpatialGraph` bundles four pieces:
#' \describe{
#'   \item{`nodes`}{a `data.table` with one row per cell. Required columns:
#'     `cell_id` (character), `x`, `y` (numeric), `cell_type` (character),
#'     `sample_id` (character). Optional: `z` (numeric, may be `NA`),
#'     `region_id` (character). Arbitrary additional columns are passed
#'     through.}
#'   \item{`edges`}{a `data.table` with one row per undirected edge. Required
#'     columns: `source`, `target` (character cell ids), `sample_id`. Optional:
#'     `distance` (numeric). Edges are stored once with `source < target`.}
#'   \item{`meta`}{a named list of build-time metadata
#'     (`build_method`, `k`, `radius`, `max_distance`, `version`, …).}
#'   \item{cache}{an internal environment used by `as_igraph()` to memoize
#'     the igraph view; not user-facing.}
#' }
#'
#' Use the constructor `SpatialGraph()` to validate and assemble these.
#'
#' @name SpatialGraph-class
#' @keywords internal
NULL

# Required schemas (additional columns are allowed and preserved).
.required_node_cols <- c("cell_id", "x", "y", "cell_type", "sample_id")
.required_edge_cols <- c("source", "target", "sample_id")

#' Construct a `SpatialGraph`
#'
#' Validates schemas, normalizes column types, ensures undirected edges are
#' canonical (`source < target`), and returns an S3 object of class
#' `SpatialGraph`. The igraph view is built lazily by `as_igraph()`.
#'
#' @param nodes a `data.frame` / `data.table` of cells (see "Details").
#' @param edges a `data.frame` / `data.table` of edges (see "Details").
#' @param meta named list of graph-construction metadata.
#' @return a `SpatialGraph`.
#' @examples
#' nodes <- data.frame(
#'   cell_id   = c("a","b","c","d"),
#'   x         = c(0,1,2,3),
#'   y         = c(0,0,1,1),
#'   cell_type = c("T","B","T","B"),
#'   sample_id = "s1",
#'   stringsAsFactors = FALSE
#' )
#' edges <- data.frame(
#'   source = c("a","b","c"), target = c("b","c","d"),
#'   distance = c(1, sqrt(2), 1), sample_id = "s1",
#'   stringsAsFactors = FALSE
#' )
#' sg <- SpatialGraph(nodes, edges, meta = list(build_method = "manual"))
#' @export
SpatialGraph <- function(nodes, edges, meta = list()) {
  if (!is.list(meta)) stop("`meta` must be a list", call. = FALSE)

  nodes <- data.table::as.data.table(nodes)
  edges <- data.table::as.data.table(edges)

  miss_n <- setdiff(.required_node_cols, names(nodes))
  if (length(miss_n)) {
    stop("nodes is missing required columns: ",
         paste(miss_n, collapse = ", "), call. = FALSE)
  }
  miss_e <- setdiff(.required_edge_cols, names(edges))
  if (length(miss_e)) {
    stop("edges is missing required columns: ",
         paste(miss_e, collapse = ", "), call. = FALSE)
  }

  # Coerce identifier columns to character — parquet readers occasionally
  # return them as factors when the file was written by a non-arrow tool.
  nodes[, cell_id := as.character(cell_id)]
  nodes[, cell_type := as.character(cell_type)]
  nodes[, sample_id := as.character(sample_id)]
  if (!"z" %in% names(nodes)) nodes[, z := NA_real_]
  edges[, source := as.character(source)]
  edges[, target := as.character(target)]
  edges[, sample_id := as.character(sample_id)]

  if (anyDuplicated(nodes$cell_id)) {
    stop("nodes$cell_id contains duplicates", call. = FALSE)
  }

  # Drop self-loops, canonicalize order so source < target.
  if (nrow(edges)) {
    self_loop <- edges$source == edges$target
    if (any(self_loop)) {
      warning(sum(self_loop), " self-loop(s) dropped", call. = FALSE)
      edges <- edges[!self_loop]
    }
    swap <- edges$source > edges$target
    if (any(swap)) {
      tmp <- edges$source[swap]
      edges[swap, source := target]
      edges[swap, target := tmp]
    }
    if (anyDuplicated(edges, by = c("source", "target"))) {
      edges <- unique(edges, by = c("source", "target"))
    }

    # Endpoint cell ids must exist in nodes.
    bad <- !(edges$source %in% nodes$cell_id) |
           !(edges$target %in% nodes$cell_id)
    if (any(bad)) {
      stop(sum(bad), " edge endpoint(s) not found in nodes$cell_id",
           call. = FALSE)
    }
  }

  data.table::setkey(nodes, cell_id)
  if (nrow(edges)) data.table::setkey(edges, source, target)

  meta_full <- utils::modifyList(
    list(build_method = "unknown", version = "0.1.0"),
    meta
  )

  structure(
    list(
      nodes  = nodes,
      edges  = edges,
      meta   = meta_full,
      .cache = new.env(parent = emptyenv())
    ),
    class = "SpatialGraph"
  )
}

#' @rdname SpatialGraph
#' @param x a `SpatialGraph`.
#' @param ... ignored.
#' @return `x` invisibly (for `print`).
#' @export
print.SpatialGraph <- function(x, ...) {
  N <- nrow(x$nodes)
  E <- nrow(x$edges)
  cat("<SpatialGraph>\n")
  cat(sprintf("  nodes: %d  |  edges: %d  |  build: %s\n",
              N, E, x$meta$build_method %||% "?"))

  samp <- sort(unique(x$nodes$sample_id))
  cat(sprintf("  samples (%d): %s\n",
              length(samp),
              paste(utils::head(samp, 5), collapse = ", ")))

  if ("region_id" %in% names(x$nodes)) {
    rg <- x$nodes[, .N, by = region_id][order(-N)]
    pretty <- paste(sprintf("%s=%d", rg$region_id, rg$N), collapse = ", ")
    cat(sprintf("  regions (%d): %s\n", nrow(rg), pretty))
  }

  ct <- x$nodes[, .N, by = cell_type][order(-N)]
  cat(sprintf("  cell types (%d, top 5): %s\n",
              nrow(ct),
              paste(sprintf("%s=%d", utils::head(ct$cell_type, 5),
                            utils::head(ct$N, 5)),
                    collapse = ", ")))

  if (!is.null(x$.cache$igraph)) {
    cat("  igraph view: cached\n")
  }
  invisible(x)
}

#' Subset a `SpatialGraph` by a node-column expression
#'
#' These wrap a non-standard-evaluation filter on `sg$nodes`. After filtering
#' nodes, edges with at least one endpoint outside the surviving set are
#' dropped. The igraph cache is invalidated.
#'
#' @param sg a `SpatialGraph`.
#' @param expr an unquoted expression evaluated against `sg$nodes` columns
#'   (e.g. `region_id == "niche_1"`, `sample_id %in% c("a", "b")`).
#' @return a new `SpatialGraph`.
#' @examples
#' nodes <- data.frame(
#'   cell_id = letters[1:5], x = 0:4, y = 0,
#'   cell_type = "T", sample_id = "s1",
#'   region_id = c("r1","r1","r2","r2","r2"),
#'   stringsAsFactors = FALSE
#' )
#' edges <- data.frame(
#'   source = c("a","b","c","d"), target = c("b","c","d","e"),
#'   sample_id = "s1", stringsAsFactors = FALSE
#' )
#' sg <- SpatialGraph(nodes, edges)
#' sg1 <- subset_region(sg, region_id == "r1")
#' nrow(sg1$nodes)  # 2
#' @export
subset_region <- function(sg, expr) {
  .subset_sg(sg, substitute(expr), parent.frame(), op = "subset_region")
}

#' @rdname subset_region
#' @export
subset_sample <- function(sg, expr) {
  .subset_sg(sg, substitute(expr), parent.frame(), op = "subset_sample")
}

.subset_sg <- function(sg, expr_quoted, env, op) {
  if (!inherits(sg, "SpatialGraph")) {
    stop("`sg` must be a SpatialGraph", call. = FALSE)
  }
  mask <- eval(expr_quoted, envir = sg$nodes, enclos = env)
  if (!is.logical(mask) || length(mask) != nrow(sg$nodes)) {
    stop("subset expression must yield a logical vector of length nrow(nodes)",
         call. = FALSE)
  }
  new_nodes <- sg$nodes[which(mask)]
  keep_ids <- new_nodes$cell_id
  new_edges <- sg$edges[source %in% keep_ids & target %in% keep_ids]

  meta <- sg$meta
  meta$build_method <- paste0(meta$build_method %||% "?", "+", op)
  SpatialGraph(new_nodes, new_edges, meta = meta)
}

#' Return (and cache) the igraph view of a `SpatialGraph`
#'
#' Builds an undirected `igraph` graph with one vertex per node and one
#' edge per row in `sg$edges`. Vertex attributes mirror `sg$nodes` columns;
#' edge attributes mirror non-endpoint columns of `sg$edges`. The graph is
#' memoized inside `sg$.cache$igraph` and reused on subsequent calls — the
#' cache is keyed on the row-counts of nodes/edges so it auto-invalidates
#' across subsets.
#'
#' @param sg a `SpatialGraph`.
#' @return an `igraph::igraph` object.
#' @examples
#' sg <- SpatialGraph(
#'   data.frame(cell_id = c("a","b"), x = 0:1, y = 0,
#'              cell_type = "T", sample_id = "s1",
#'              stringsAsFactors = FALSE),
#'   data.frame(source = "a", target = "b", sample_id = "s1",
#'              stringsAsFactors = FALSE)
#' )
#' g <- as_igraph(sg)
#' igraph::vcount(g)
#' @export
as_igraph <- function(sg) {
  if (!inherits(sg, "SpatialGraph")) {
    stop("`sg` must be a SpatialGraph", call. = FALSE)
  }
  key <- c(nrow(sg$nodes), nrow(sg$edges))
  cached <- sg$.cache$igraph
  if (!is.null(cached) && identical(sg$.cache$key, key)) {
    return(cached)
  }
  vertex_attr <- as.data.frame(sg$nodes)
  edge_attr <- if (nrow(sg$edges)) {
    as.data.frame(sg$edges)
  } else {
    data.frame(source = character(), target = character(),
               stringsAsFactors = FALSE)
  }
  g <- igraph::graph_from_data_frame(
    d = edge_attr,
    directed = FALSE,
    vertices = vertex_attr
  )
  sg$.cache$igraph <- g
  sg$.cache$key <- key
  g
}
