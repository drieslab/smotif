#' Build a `SpatialGraph` from coordinates
#'
#' Build a `SpatialGraph` over a set of cells given their spatial coordinates,
#' cell types, and sample id, using one of three connectivity rules:
#' \itemize{
#'   \item `"knn"` (default): each cell connects to its `k` nearest neighbors.
#'     Edges are deduplicated and stored undirected, so the resulting degree
#'     can exceed `k` (a cell may be in another's k-NN without the converse).
#'   \item `"delaunay"`: edges from the 2D Delaunay triangulation (xy plane).
#'   \item `"radius"`: edges between every pair within `radius` of each other.
#' }
#'
#' Implementations rely on lightweight optional dependencies (`FNN`, `deldir`,
#' `dbscan`); install only the one(s) you need. An optional `max_distance`
#' cap is applied after edge construction.
#'
#' @param coords numeric matrix or `data.frame` with columns `x, y` and
#'   optionally `z` (z is preserved on nodes but ignored by all
#'   connectivity rules in v0.1).
#' @param cell_types character vector aligned with rows of `coords`.
#' @param sample_id character; recycled to length `nrow(coords)` if scalar.
#' @param region_id optional character vector aligned with rows of `coords`.
#' @param method one of `"knn"`, `"delaunay"`, `"radius"`.
#' @param k integer; neighbors when `method = "knn"`.
#' @param radius numeric; cutoff when `method = "radius"`.
#' @param max_distance optional numeric cap applied after edge construction.
#' @param cell_id optional character vector of cell ids; if missing, ids
#'   `cell_<i>` are generated.
#' @return a `SpatialGraph`.
#' @examples
#' set.seed(1)
#' coords <- cbind(x = runif(20), y = runif(20))
#' sg <- build_spatial_graph(
#'   coords, cell_types = sample(c("T","B"), 20, TRUE),
#'   sample_id = "s1", method = "knn", k = 3
#' )
#' nrow(sg$edges)
#' @export
build_spatial_graph <- function(coords,
                                cell_types,
                                sample_id,
                                region_id = NULL,
                                method = c("knn", "delaunay", "radius"),
                                k = 6L,
                                radius = NULL,
                                max_distance = NULL,
                                cell_id = NULL) {
  method <- match.arg(method)
  coords <- .coerce_coords(coords)
  N <- nrow(coords)
  if (length(cell_types) != N) {
    stop("length(cell_types) must equal nrow(coords)", call. = FALSE)
  }
  if (length(sample_id) == 1L) sample_id <- rep(sample_id, N)
  if (length(sample_id) != N) {
    stop("sample_id must be length 1 or nrow(coords)", call. = FALSE)
  }
  if (is.null(cell_id)) cell_id <- paste0("cell_", seq_len(N))
  if (length(cell_id) != N || anyDuplicated(cell_id)) {
    stop("cell_id must be length nrow(coords) with unique values",
         call. = FALSE)
  }
  cell_id <- as.character(cell_id)

  edges_idx <- switch(
    method,
    knn      = .edges_knn(coords[, c("x", "y"), drop = FALSE], k),
    radius   = {
      if (is.null(radius)) {
        stop("method = 'radius' requires a non-null `radius`", call. = FALSE)
      }
      .edges_radius(coords[, c("x", "y"), drop = FALSE], radius)
    },
    delaunay = .edges_delaunay(coords[, c("x", "y"), drop = FALSE])
  )

  edges <- data.table::data.table(
    source    = cell_id[edges_idx$source],
    target    = cell_id[edges_idx$target],
    distance  = edges_idx$distance,
    sample_id = sample_id[edges_idx$source]
  )
  if (!is.null(max_distance)) {
    edges <- edges[distance <= max_distance]
  }

  nodes <- data.table::data.table(
    cell_id   = cell_id,
    x         = coords[, "x"],
    y         = coords[, "y"],
    z         = if ("z" %in% colnames(coords)) coords[, "z"] else NA_real_,
    cell_type = as.character(cell_types),
    sample_id = sample_id
  )
  if (!is.null(region_id)) {
    if (length(region_id) != N) {
      stop("region_id must be length nrow(coords)", call. = FALSE)
    }
    nodes[, region_id := as.character(region_id)]
  }

  meta <- list(
    build_method = method,
    k            = if (method == "knn") k else NA_integer_,
    radius       = if (method == "radius") radius else NA_real_,
    max_distance = max_distance %||% NA_real_,
    n_nodes      = N
  )
  SpatialGraph(nodes, edges, meta = meta)
}

# ---- internal: coordinate coercion -----------------------------------------

.coerce_coords <- function(coords) {
  if (is.data.frame(coords) || data.table::is.data.table(coords)) {
    cols <- intersect(c("x", "y", "z"), names(coords))
    if (!all(c("x", "y") %in% cols)) {
      stop("coords must have columns 'x' and 'y'", call. = FALSE)
    }
    coords <- as.matrix(coords[, cols, drop = FALSE])
  } else {
    coords <- as.matrix(coords)
    if (is.null(colnames(coords))) {
      if (ncol(coords) < 2L) {
        stop("coords matrix must have at least 2 columns", call. = FALSE)
      }
      colnames(coords) <- c("x", "y", if (ncol(coords) >= 3L) "z")[
        seq_len(min(ncol(coords), 3L))
      ]
    }
    if (!all(c("x", "y") %in% colnames(coords))) {
      stop("coords must have columns 'x' and 'y'", call. = FALSE)
    }
  }
  storage.mode(coords) <- "double"
  coords
}

# ---- internal: connectivity rules ------------------------------------------

# Each .edges_* helper returns a list with integer source/target column
# vectors (1-based row indices into coords) and a numeric `distance` vector.
# Endpoints are canonicalized (source < target) and deduplicated here so
# the SpatialGraph constructor does not have to redo that work.

.require_suggest <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "Package '%s' is required for %s. Install it with install.packages('%s').",
      pkg, why, pkg), call. = FALSE)
  }
}

.canonicalize_edges <- function(s, t, d) {
  swap <- s > t
  s2 <- ifelse(swap, t, s)
  t2 <- ifelse(swap, s, t)
  dt <- data.table::data.table(s = s2, t = t2, distance = d)
  dt <- unique(dt, by = c("s", "t"))
  list(source = dt$s, target = dt$t, distance = dt$distance)
}

.edges_knn <- function(xy, k) {
  .require_suggest("FNN", "method = 'knn'")
  if (nrow(xy) <= k) {
    stop("nrow(coords) must be greater than k for method = 'knn'",
         call. = FALSE)
  }
  nn <- FNN::get.knn(xy, k = k)
  N  <- nrow(xy)
  s  <- rep(seq_len(N), each = k)
  t  <- as.integer(t(nn$nn.index))
  d  <- as.numeric(t(nn$nn.dist))
  .canonicalize_edges(s, t, d)
}

.edges_radius <- function(xy, radius) {
  .require_suggest("dbscan", "method = 'radius'")
  fr <- dbscan::frNN(xy, eps = radius)
  # frNN returns per-row id and dist lists. Flatten to source/target.
  lens <- lengths(fr$id)
  if (sum(lens) == 0L) {
    return(list(source = integer(), target = integer(), distance = numeric()))
  }
  s <- rep(seq_along(fr$id), times = lens)
  t <- unlist(fr$id, use.names = FALSE)
  d <- unlist(fr$dist, use.names = FALSE)
  .canonicalize_edges(s, t, d)
}

.edges_delaunay <- function(xy) {
  .require_suggest("deldir", "method = 'delaunay'")
  d <- deldir::deldir(xy[, "x"], xy[, "y"], suppressMsge = TRUE)
  s <- d$delsgs$ind1
  t <- d$delsgs$ind2
  dist <- sqrt((xy[s, "x"] - xy[t, "x"])^2 + (xy[s, "y"] - xy[t, "y"])^2)
  .canonicalize_edges(s, t, dist)
}
