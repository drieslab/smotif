#' Read a `SpatialGraph` from parquet
#'
#' Reads `nodes.parquet` and `edges.parquet`, validates the schemas, and
#' returns a `SpatialGraph`. If a sidecar `<nodes_path>.meta.json` exists,
#' it is loaded into `sg$meta`; otherwise `meta` records the source paths
#' so a round-trip is reproducible.
#'
#' @param nodes_path path to a `nodes.parquet` file.
#' @param edges_path path to an `edges.parquet` file.
#' @return a `SpatialGraph`.
#' @examples
#' \dontrun{
#'   sg <- read_spatial_graph("nodes.parquet", "edges.parquet")
#' }
#' @export
read_spatial_graph <- function(nodes_path, edges_path) {
  if (!file.exists(nodes_path)) stop("nodes_path not found: ", nodes_path,
                                     call. = FALSE)
  if (!file.exists(edges_path)) stop("edges_path not found: ", edges_path,
                                     call. = FALSE)
  nodes <- arrow::read_parquet(nodes_path)
  edges <- arrow::read_parquet(edges_path)
  meta <- list(
    build_method = "read_spatial_graph",
    nodes_path   = normalizePath(nodes_path, mustWork = TRUE),
    edges_path   = normalizePath(edges_path, mustWork = TRUE)
  )
  sidecar <- paste0(nodes_path, ".meta.json")
  if (file.exists(sidecar) && requireNamespace("jsonlite", quietly = TRUE)) {
    extra <- tryCatch(
      jsonlite::fromJSON(sidecar, simplifyVector = TRUE),
      error = function(e) NULL
    )
    if (is.list(extra)) meta <- utils::modifyList(meta, extra)
  }
  SpatialGraph(nodes, edges, meta = meta)
}

#' Write a `SpatialGraph` to parquet
#'
#' Writes `sg$nodes` and `sg$edges` as separate parquet files. If
#' `jsonlite` is available, also writes `<nodes_path>.meta.json` so that
#' `meta` survives a round trip.
#'
#' @param sg a `SpatialGraph`.
#' @param nodes_path output path for `nodes.parquet`.
#' @param edges_path output path for `edges.parquet`.
#' @return `sg`, invisibly.
#' @examples
#' \dontrun{
#'   write_spatial_graph(sg, "nodes.parquet", "edges.parquet")
#' }
#' @export
write_spatial_graph <- function(sg, nodes_path, edges_path) {
  if (!inherits(sg, "SpatialGraph")) {
    stop("`sg` must be a SpatialGraph", call. = FALSE)
  }
  arrow::write_parquet(sg$nodes, nodes_path)
  arrow::write_parquet(sg$edges, edges_path)
  if (length(sg$meta) && requireNamespace("jsonlite", quietly = TRUE)) {
    sidecar <- paste0(nodes_path, ".meta.json")
    # Strip any environments / non-JSON-friendly fields out of meta first.
    safe_meta <- sg$meta[vapply(sg$meta, function(v) {
      is.atomic(v) || is.list(v)
    }, logical(1))]
    jsonlite::write_json(safe_meta, sidecar,
                         pretty = TRUE, auto_unbox = TRUE)
  }
  invisible(sg)
}
