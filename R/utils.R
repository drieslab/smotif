#' @importFrom data.table := .N .SD data.table is.data.table
#'   setDT setkey setnames as.data.table rbindlist
#' @keywords internal
"_PACKAGE"

# Declare names referenced by data.table NSE so that `R CMD check` does
# not flag them as undefined globals. Each name appears in a `dt[i = ...,
# j = ...]` expression somewhere in the package.
utils::globalVariables(c(
  ".N", ".SD", ".I",
  "cell_id", "cell_type", "region_id", "sample_id",
  "source", "target", "distance",
  "leiden_0.4", "leiden_0.8",
  "z", "N", "cell_ID",
  # motif backends
  "ca", "cb", "cc", "c1", "c2", "c3",
  "color", "color_tuple", "color_lo", "color_hi",
  "color_min", "color_mid", "color_max",
  "deg", "instance_id", "motif_id", "orbit_id",
  "x_a", "x_c", "end_lo", "end_hi", "ac_key", "count",
  # build_graph internals
  "s", "t",
  # catalog columns referenced via NSE in .cap_instances
  "size", "canonical_iso_class",
  ".",
  # wedge / triangle vertex columns
  "v1", "v2", "v3", "b",
  # gene-association internals
  "in_cells", "n_in", "n_out", "p", "p_adj_BH",
  "log2FC", "z", "z_combined", "p_combined",
  "all_in",
  # compare_motifs internals
  "total", "chi2", "df",
  # size-4 backend internals
  "w1", "w2", "w3", "w4", "d1", "d2", "d3", "d4",
  "ds1", "ds2", "ds3", "ds4",
  "e12", "e13", "e14", "e23", "e24", "e34",
  "n_edges", "x_a", "x_c",
  "x", "d", "c4"
))

#' Null-coalescing operator
#'
#' Returns `b` when `a` is `NULL` or zero-length, else `a`.
#' @param a,b values
#' @return `a` if non-null and non-empty, otherwise `b`
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

#' Stable md5 hash of an arbitrary R object
#'
#' Used for cache keys (e.g., LLM marker-input fingerprinting). Computes
#' md5 over `serialize()` so we don't pull in another dependency.
#' @param x any R object
#' @return character scalar, hex md5 digest
#' @keywords internal
#' @noRd
.smotif_hash <- function(x) {
  tf <- tempfile()
  on.exit(unlink(tf), add = TRUE)
  con <- file(tf, open = "wb")
  serialize(x, con)
  close(con)
  unname(tools::md5sum(tf))
}

#' Tiny timestamped console logger
#'
#' Used by `inst/scripts/prepare_workshop_xenium.R` and by long-running
#' internal routines to make progress visible without pulling in a
#' logging dependency.
#' @param msg character scalar
#' @keywords internal
#' @noRd
.step_log <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  invisible(NULL)
}

# Internal: touch one Matrix symbol so `R CMD check` doesn't flag the
# Imports entry as unused while sparse-matrix code lands in step 6.
# Removed once `genes.R` references Matrix directly.
.imports_touch <- function() {
  Matrix::sparseMatrix(i = integer(), j = integer(),
                       x = numeric(), dims = c(0L, 0L))
}
