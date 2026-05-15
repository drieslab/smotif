#' Test motif enrichment under a spatially-aware null
#'
#' Compares observed per-motif counts against a permutation null
#' distribution and returns a `data.table` of test statistics. Three null
#' models are supported (see [null_models][smotif::test_motif_enrichment]):
#'
#' \itemize{
#'   \item `"label_permutation"`: shuffle `cell_type` across all nodes;
#'     graph topology stays fixed.
#'   \item `"region_stratified"`: shuffle `cell_type` within each
#'     `region_id`; preserves per-region cell-type composition. Requires
#'     `sg$nodes$region_id`.
#'   \item `"geometric_jitter"`: perturb coordinates by Gaussian noise of
#'     standard deviation `null_args$sd`, rebuild the graph using the
#'     original `build_spatial_graph` parameters, recount. Requires
#'     `sg$meta$build_method` to be one of `"knn"`, `"delaunay"`,
#'     `"radius"`.
#' }
#'
#' Empirical p-values use the unbiased "+1" estimator
#' \eqn{p = (1 + \sum I[null_i \ge observed]) / (1 + n_{perm})} so the
#' result is in `(0, 1]` for finite `n_perm`. The fold change is
#' `observed / max(mean_null, .Machine$double.eps)`. Adjustment for
#' multiple testing is BH on `p_emp`.
#'
#' @param motifs a `Motifs` object from [find_motifs()].
#' @param sg the `SpatialGraph` `motifs` was computed from.
#' @param null one of `"label_permutation"`, `"region_stratified"`,
#'   `"geometric_jitter"`.
#' @param n_perm integer; number of permutations.
#' @param seed integer; RNG seed for reproducibility.
#' @param parallel logical; future hook — v0.1 ignores and runs serially.
#' @param null_args optional list of arguments specific to the chosen
#'   null. Currently used only by `"geometric_jitter"` (`sd`).
#' @return a `data.table` with columns
#'   `motif_id, observed, mean_null, sd_null, z, p_emp, p_adj_BH, fold`.
#' @examples
#' set.seed(1)
#' coords <- cbind(x = runif(40), y = runif(40))
#' sg <- build_spatial_graph(
#'   coords, cell_types = sample(c("T","B"), 40, TRUE),
#'   sample_id = "s1", method = "knn", k = 4
#' )
#' m <- find_motifs(sg, size = 3L)
#' test_motif_enrichment(m, sg, null = "label_permutation",
#'                       n_perm = 50L, seed = 1L)
#' @export
test_motif_enrichment <- function(motifs,
                                  sg,
                                  null = c("label_permutation",
                                           "region_stratified",
                                           "geometric_jitter"),
                                  n_perm = 1000L,
                                  seed = 1L,
                                  parallel = FALSE,
                                  null_args = list()) {
  if (!inherits(motifs, "Motifs")) {
    stop("`motifs` must be a Motifs object", call. = FALSE)
  }
  if (!inherits(sg, "SpatialGraph")) {
    stop("`sg` must be a SpatialGraph", call. = FALSE)
  }
  null <- match.arg(null)
  n_perm <- as.integer(n_perm)
  if (length(n_perm) != 1L || is.na(n_perm) || n_perm < 1L) {
    stop("n_perm must be a positive integer", call. = FALSE)
  }
  if (isTRUE(parallel)) {
    # TODO: wire up parallel::mclapply / future.apply backend in v0.2.
    warning("parallel = TRUE is not implemented in v0.1; running serially",
            call. = FALSE)
  }

  observed_dt <- motifs$catalog[, .(motif_id, observed = count)]
  obs_ids <- observed_dt$motif_id

  if (null %in% c("label_permutation", "region_stratified")) {
    if (length(motifs$instance_meta) == 0L) {
      stop(
        "label permutation null needs `motifs$instance_meta`. ",
        "It is populated by find_motifs() but missing on hand-built ",
        "Motifs objects.", call. = FALSE
      )
    }
  }

  set.seed(seed)
  draw <- switch(null,
    label_permutation = function() .draw_label_perm(motifs, sg, obs_ids,
                                                    stratified = FALSE,
                                                    null_args),
    region_stratified = function() .draw_label_perm(motifs, sg, obs_ids,
                                                    stratified = TRUE,
                                                    null_args),
    geometric_jitter  = function() .draw_jitter(motifs, sg, obs_ids,
                                                null_args)
  )

  # null_mat: rows = motifs, cols = permutations.
  null_mat <- matrix(0L, nrow = length(obs_ids), ncol = n_perm,
                     dimnames = list(obs_ids, NULL))
  for (i in seq_len(n_perm)) null_mat[, i] <- draw()

  observed <- observed_dt$observed
  names(observed) <- obs_ids
  mean_null <- rowMeans(null_mat)
  sd_null   <- apply(null_mat, 1L, stats::sd)
  z <- (observed - mean_null) /
    pmax(sd_null, sqrt(.Machine$double.eps))
  ge <- rowSums(null_mat >= observed)
  p_emp <- (1 + ge) / (1 + n_perm)
  p_adj <- stats::p.adjust(p_emp, method = "BH")
  fold  <- observed / pmax(mean_null, .Machine$double.eps)

  out <- data.table::data.table(
    motif_id  = obs_ids,
    observed  = observed,
    mean_null = mean_null,
    sd_null   = sd_null,
    z         = z,
    p_emp     = p_emp,
    p_adj_BH  = p_adj,
    fold      = fold
  )
  data.table::setorder(out, p_emp, -z)
  out[]
}
