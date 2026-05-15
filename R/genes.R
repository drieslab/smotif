#' Associate genes (or gene-set scores) with motifs or motif orbits
#'
#' For each motif (or motif × orbit), test whether genes are differentially
#' expressed between cells participating in that motif/orbit ("in") and
#' the rest of the dataset ("out"). Two methods are supported:
#'
#' \itemize{
#'   \item `"wilcoxon"`: rank-sum test per gene with tie correction. The
#'     full per-cell rank matrix is computed once and reused across motifs.
#'     Genes are filtered by `min_in_expressed` (number of "in" cells with
#'     non-zero expression).
#'   \item `"module_score"`: per-cell scores are computed as the mean of
#'     row-standardized (z-scored) expression over each gene set, then
#'     compared between motif and background using the same rank-sum
#'     framework.
#' }
#'
#' Two pooling modes are supported:
#' \itemize{
#'   \item `pooling = "pooled"` (default): all cells across all samples
#'     contribute to a single in/out test per motif. Fastest; appropriate
#'     for single-sample data or when sample-level heterogeneity is not
#'     of interest.
#'   \item `pooling = "per_sample"`: each sample is tested independently
#'     and per-sample statistics are combined via a Stouffer-weighted Z
#'     meta-analysis (weights `= sqrt(n_in)` per sample). The output
#'     reports `n_samples_tested`, `mean_log2FC` and `sd_log2FC` so
#'     between-sample heterogeneity is visible: a gene with a large
#'     `z_combined` and a tiny `sd_log2FC` reproduces across samples; a
#'     gene with a large `z_combined` and large `sd_log2FC` is being
#'     driven by one sample.
#' }
#'
#' @param motifs a `Motifs` object from [find_motifs()].
#' @param expression matrix or `Matrix::dgCMatrix`, genes × cells, with
#'   row names (genes) and column names (cell ids matching `sg`).
#' @param sg the `SpatialGraph` `motifs` was computed from.
#' @param by `"motif"` (one test per motif) or `"orbit"` (one test per
#'   motif × orbit_id combination).
#' @param method `"wilcoxon"` or `"module_score"`.
#' @param modules named list of gene sets; required when
#'   `method = "module_score"`.
#' @param min_cells motifs/orbits with fewer than this many unique cells
#'   in the "in" set are skipped (per-sample under `pooling = "per_sample"`).
#' @param min_in_expressed for `method = "wilcoxon"`, only test genes
#'   expressed in at least this many "in" cells.
#' @param pooling `"pooled"` (default) or `"per_sample"`. See "Details".
#' @return a `data.table`. Columns under `pooling = "pooled"`:
#'   `motif_id, orbit_id, gene, n_in, n_out, n_in_expressed,
#'   mean_in, mean_out, log2FC, U, z, p, p_adj_BH`. Under
#'   `pooling = "per_sample"`: `motif_id, orbit_id, gene,
#'   n_samples_tested, n_in_total, n_out_total, mean_log2FC, sd_log2FC,
#'   z_combined, p_combined, p_adj_BH`.
#' @examples
#' set.seed(1)
#' coords <- cbind(x = runif(60), y = runif(60))
#' types  <- sample(c("T", "B"), 60, TRUE)
#' sg <- build_spatial_graph(coords, types, "s1", method = "knn", k = 4L)
#' m  <- find_motifs(sg, size = 3L)
#' expr <- Matrix::Matrix(
#'   matrix(stats::rpois(60 * 50, 1), nrow = 50),
#'   sparse = TRUE
#' )
#' rownames(expr) <- paste0("g", 1:50)
#' colnames(expr) <- sg$nodes$cell_id
#' associate_genes(m, expr, sg, min_cells = 8L)[1:5]
#' @export
associate_genes <- function(motifs,
                            expression,
                            sg,
                            by = c("motif", "orbit"),
                            method = c("wilcoxon", "module_score"),
                            modules = NULL,
                            min_cells = 10L,
                            min_in_expressed = 3L,
                            pooling = c("pooled", "per_sample")) {
  by <- match.arg(by)
  method <- match.arg(method)
  pooling <- match.arg(pooling)
  if (!inherits(motifs, "Motifs")) {
    stop("`motifs` must be a Motifs object", call. = FALSE)
  }
  if (!inherits(sg, "SpatialGraph")) {
    stop("`sg` must be a SpatialGraph", call. = FALSE)
  }
  if (!is.matrix(expression) && !inherits(expression, "Matrix")) {
    stop("`expression` must be a matrix or Matrix (sparse)", call. = FALSE)
  }
  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("`expression` must have row names (genes) and column names (cell ids)",
         call. = FALSE)
  }
  if (identical(method, "module_score") &&
      (!is.list(modules) || !length(modules))) {
    stop('method = "module_score" requires a non-empty `modules` list',
         call. = FALSE)
  }

  if (identical(pooling, "per_sample")) {
    return(.assoc_per_sample(motifs, expression, sg, by, method,
                             modules, min_cells, min_in_expressed))
  }

  groups <- .group_in_cells(motifs, expression, by, min_cells)
  if (!nrow(groups)) {
    return(.empty_assoc_dt(method))
  }

  if (identical(method, "wilcoxon")) {
    .assoc_wilcoxon(groups, expression, by, min_in_expressed)
  } else {
    .assoc_module_score(groups, expression, modules, by)
  }
}

# ---- internal: group definition ------------------------------------------

# For each (motif_id [+ orbit_id]) group, return its unique "in" cell ids
# restricted to cells that appear in the expression columns. Skip groups
# with fewer than `min_cells` unique cells.
.group_in_cells <- function(motifs, expression, by, min_cells) {
  cell_ids <- colnames(expression)
  inc <- motifs$incidence[cell_id %in% cell_ids]
  if (!nrow(inc)) {
    return(data.table::data.table(
      motif_id = character(), orbit_id = integer(),
      in_cells = list(), n_in = integer()
    ))
  }
  group_cols <- if (identical(by, "orbit")) c("motif_id", "orbit_id") else "motif_id"
  groups <- inc[, .(in_cells = list(unique(cell_id))), by = group_cols]
  groups[, n_in := vapply(in_cells, length, integer(1L))]
  groups <- groups[n_in >= min_cells]
  if (!"orbit_id" %in% names(groups)) groups[, orbit_id := NA_integer_]
  groups
}

# ---- internal: shared rank-sum core --------------------------------------

# Per-row rank with ties resolved by averaging. Returns a (rows × ncol)
# numeric matrix the same shape as the input. For sparse input we densify
# first; v0.1 caveat: for very wide expression panels (>~ 20k genes by
# >~ 50k cells) this can blow memory — densification is fine for the
# workshop dataset (397 × 7586 ≈ 24 MB).
.rank_per_row <- function(x_dense) {
  out <- matrix(0.0, nrow = nrow(x_dense), ncol = ncol(x_dense),
                dimnames = dimnames(x_dense))
  for (g in seq_len(nrow(x_dense))) {
    out[g, ] <- rank(x_dense[g, ])
  }
  out
}

# Per-row tie-correction term sum(t^3 - t) over groups of tied values.
# Used to deflate the rank-sum variance when there are many ties (the
# usual situation for sparse expression: huge tied group at zero).
.tie_term_per_row <- function(x_dense) {
  vapply(seq_len(nrow(x_dense)), function(g) {
    tab <- tabulate(match(x_dense[g, ], unique(x_dense[g, ])))
    tg <- tab[tab > 1L]
    if (!length(tg)) 0 else sum(tg^3 - tg)
  }, numeric(1L))
}

# Mann-Whitney U / rank-sum two-sided test, vectorized over rows.
# Inputs:
#   rank_mat   numeric matrix, ranks of each gene over all cells
#   row_sums   numeric vector, total ranked-row sum per row (i.e. N(N+1)/2)
#   tie_terms  numeric vector, per-row tie correction
#   in_idx     integer vector, columns belonging to the "in" set
#   N          total number of cells
# Returns a list of equal-length numeric vectors: U, z, p.
.rank_sum_test <- function(rank_mat, row_sums, tie_terms, in_idx, N) {
  n_in <- length(in_idx)
  n_out <- N - n_in
  R1 <- rowSums(rank_mat[, in_idx, drop = FALSE])
  U  <- R1 - n_in * (n_in + 1) / 2
  mu <- n_in * n_out / 2
  sigma2 <- (n_in * n_out / 12) * (
    (N + 1) - tie_terms / (N * (N - 1))
  )
  sigma <- sqrt(pmax(sigma2, .Machine$double.eps))
  z <- (U - mu) / sigma
  p <- 2 * (1 - stats::pnorm(abs(z)))
  list(U = unname(U), z = unname(z), p = unname(p))
}

# ---- internal: Wilcoxon over genes ---------------------------------------

.assoc_wilcoxon <- function(groups, expression, by, min_in_expressed) {
  cell_ids <- colnames(expression)
  N <- length(cell_ids)
  G <- nrow(expression)

  # Densify once. See .rank_per_row for the memory caveat.
  expr_dense <- if (is.matrix(expression)) {
    expression
  } else {
    as.matrix(expression)
  }
  rank_mat  <- .rank_per_row(expr_dense)
  tie_terms <- .tie_term_per_row(expr_dense)
  row_sums  <- rowSums(expr_dense)

  results <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    grp <- groups[i]
    in_cells <- grp$in_cells[[1L]]
    n_in <- grp$n_in
    n_out <- N - n_in
    in_idx <- match(in_cells, cell_ids)
    in_idx <- in_idx[!is.na(in_idx)]
    if (length(in_idx) < grp$n_in) n_in <- length(in_idx)
    if (n_in < 1L || n_out < 1L) next

    test <- .rank_sum_test(rank_mat, row_sums, tie_terms, in_idx, N)

    sub_in <- expr_dense[, in_idx, drop = FALSE]
    sum_in <- rowSums(sub_in)
    mean_in  <- sum_in / n_in
    mean_out <- (row_sums - sum_in) / n_out
    log2fc <- log2((mean_in + 1) / (mean_out + 1))
    n_in_expr <- rowSums(sub_in > 0)

    keep <- n_in_expr >= min_in_expressed
    if (!any(keep)) next

    results[[i]] <- data.table::data.table(
      motif_id       = grp$motif_id,
      orbit_id       = grp$orbit_id,
      gene           = rownames(expression)[keep],
      n_in           = n_in,
      n_out          = n_out,
      n_in_expressed = n_in_expr[keep],
      mean_in        = mean_in[keep],
      mean_out       = mean_out[keep],
      log2FC         = log2fc[keep],
      U              = test$U[keep],
      z              = test$z[keep],
      p              = test$p[keep]
    )
  }

  out <- data.table::rbindlist(results)
  if (nrow(out)) {
    out[, p_adj_BH := stats::p.adjust(p, method = "BH"), by = motif_id]
    data.table::setorder(out, motif_id, p)
  } else {
    out <- .empty_assoc_dt("wilcoxon")
  }
  out
}

# ---- internal: module score ---------------------------------------------

.assoc_module_score <- function(groups, expression, modules, by) {
  cell_ids <- colnames(expression)
  N <- length(cell_ids)

  expr_dense <- if (is.matrix(expression)) expression else as.matrix(expression)

  # Row-wise z-score (each gene standardized across cells) using a numeric
  # floor on sd to avoid divide-by-zero for genes with no expression.
  mu <- rowMeans(expr_dense)
  msq <- rowMeans(expr_dense * expr_dense)
  sd_g <- sqrt(pmax(msq - mu^2, 0))
  sd_g <- pmax(sd_g, sqrt(.Machine$double.eps))
  z_expr <- (expr_dense - mu) / sd_g

  # Module × cell score matrix.
  mod_names <- names(modules) %||% paste0("module_", seq_along(modules))
  M <- length(modules)
  mod_mat <- matrix(NA_real_, nrow = M, ncol = N,
                    dimnames = list(mod_names, cell_ids))
  for (k in seq_along(modules)) {
    genes <- intersect(modules[[k]], rownames(expr_dense))
    if (!length(genes)) next
    mod_mat[k, ] <- colMeans(z_expr[genes, , drop = FALSE])
  }
  keep_mod <- !apply(is.na(mod_mat), 1L, all)
  mod_mat <- mod_mat[keep_mod, , drop = FALSE]
  if (!nrow(mod_mat)) return(.empty_assoc_dt("module_score"))

  rank_mat  <- .rank_per_row(mod_mat)
  tie_terms <- .tie_term_per_row(mod_mat)
  row_sums  <- rowSums(mod_mat)

  results <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    grp <- groups[i]
    in_cells <- grp$in_cells[[1L]]
    in_idx <- match(in_cells, cell_ids)
    in_idx <- in_idx[!is.na(in_idx)]
    n_in  <- length(in_idx); n_out <- N - n_in
    if (n_in < 1L || n_out < 1L) next

    test <- .rank_sum_test(rank_mat, row_sums, tie_terms, in_idx, N)

    sub_in <- mod_mat[, in_idx, drop = FALSE]
    mean_in  <- rowMeans(sub_in)
    mean_out <- (row_sums - rowSums(sub_in)) / n_out
    # log2FC on a z-score is undefined; report difference instead.
    diff_score <- mean_in - mean_out

    results[[i]] <- data.table::data.table(
      motif_id  = grp$motif_id,
      orbit_id  = grp$orbit_id,
      gene      = rownames(mod_mat),
      n_in      = n_in,
      n_out     = n_out,
      n_in_expressed = n_in,
      mean_in   = mean_in,
      mean_out  = mean_out,
      log2FC    = diff_score,
      U         = test$U,
      z         = test$z,
      p         = test$p
    )
  }
  out <- data.table::rbindlist(results)
  if (nrow(out)) {
    out[, p_adj_BH := stats::p.adjust(p, method = "BH"), by = motif_id]
    data.table::setorder(out, motif_id, p)
  } else {
    out <- .empty_assoc_dt("module_score")
  }
  out
}

# ---- internal: per-sample meta-analysis ----------------------------------

# Run associate_genes(pooling = "pooled") on each sample independently
# and combine via Stouffer-weighted Z. Uses sqrt(n_in) per sample as the
# weight, the standard choice when per-sample tests have unit variance
# under H0 and effective sample size differs.
.assoc_per_sample <- function(motifs, expression, sg, by, method,
                              modules, min_cells, min_in_expressed) {
  samples <- sort(unique(sg$nodes$sample_id))
  if (length(samples) < 2L) {
    warning(
      "only one sample in sg; pooling = 'per_sample' falls back to 'pooled'",
      call. = FALSE
    )
    return(associate_genes(motifs, expression, sg,
                           by = by, method = method,
                           modules = modules, min_cells = min_cells,
                           min_in_expressed = min_in_expressed,
                           pooling = "pooled"))
  }

  per_sample <- vector("list", length(samples))
  for (i in seq_along(samples)) {
    sid <- samples[i]
    cells_s <- sg$nodes[sample_id == sid, cell_id]
    keep_cols <- intersect(colnames(expression), cells_s)
    if (!length(keep_cols)) next
    expr_s <- expression[, keep_cols, drop = FALSE]
    motifs_s <- .subset_motifs_for_sample(motifs, cells_s)
    if (!nrow(motifs_s$catalog)) next
    sg_s <- subset_sample(sg, sample_id == sid)
    res_s <- associate_genes(
      motifs_s, expr_s, sg_s,
      by = by, method = method, modules = modules,
      min_cells = min_cells, min_in_expressed = min_in_expressed,
      pooling = "pooled"
    )
    if (nrow(res_s)) {
      res_s[, sample_id := sid]
      per_sample[[i]] <- res_s
    }
  }

  long <- data.table::rbindlist(per_sample)
  if (!nrow(long)) return(.empty_assoc_per_sample_dt())
  .combine_meta(long)
}

# Restrict a Motifs object to instances whose cells all live in
# `cells_in_sample`. Edges are by construction within-sample (the prep
# script and build_spatial_graph never join across samples), so checking
# the first cell of each instance suffices; we still fall back to the
# all-cells check defensively.
.subset_motifs_for_sample <- function(motifs, cells_in_sample) {
  if (!nrow(motifs$incidence)) return(motifs)
  inst_membership <- motifs$incidence[
    , .(all_in = all(cell_id %in% cells_in_sample)),
    by = instance_id
  ]
  keep_inst <- inst_membership[all_in == TRUE, instance_id]
  inc_sub <- motifs$incidence[instance_id %in% keep_inst]
  if (!nrow(inc_sub)) {
    return(Motifs(motifs$catalog[0L], motifs$incidence[0L],
                  list(), motifs$meta))
  }
  cat_sub <- inc_sub[, .(count = data.table::uniqueN(instance_id)),
                     by = motif_id]
  cat_sub <- merge(
    cat_sub,
    unique(motifs$catalog[, .(motif_id, size, canonical_iso_class,
                              color_tuple)]),
    by = "motif_id"
  )
  data.table::setcolorder(cat_sub,
    c("motif_id", "size", "canonical_iso_class", "color_tuple", "count"))
  inst_meta_sub <- list()
  for (cls in names(motifs$instance_meta)) {
    tab <- motifs$instance_meta[[cls]]
    sub <- tab[instance_id %in% keep_inst]
    if (nrow(sub)) inst_meta_sub[[cls]] <- sub
  }
  Motifs(cat_sub, inc_sub, inst_meta_sub, motifs$meta)
}

# Stouffer-weighted Z meta-analysis over per-sample rows. Each row's
# contribution is sqrt(n_in_i) * z_i / sqrt(sum n_in_i), so larger
# samples (where the per-sample test has more power) carry more weight.
# `mean_log2FC` and `sd_log2FC` track between-sample heterogeneity:
# small sd ⇒ reproducible effect, large sd ⇒ one sample driving the
# combined statistic.
.combine_meta <- function(long) {
  group_cols <- c("motif_id", "orbit_id", "gene")
  combined <- long[, .(
    n_samples_tested = .N,
    n_in_total       = sum(n_in),
    n_out_total      = sum(n_out),
    mean_log2FC      = mean(log2FC, na.rm = TRUE),
    sd_log2FC        = stats::sd(log2FC, na.rm = TRUE),
    z_combined       = sum(sqrt(n_in) * z) / sqrt(sum(n_in))
  ), by = group_cols]
  combined[, p_combined := 2 * stats::pnorm(-abs(z_combined))]
  combined[, p_adj_BH   := stats::p.adjust(p_combined, method = "BH"),
           by = motif_id]
  data.table::setorder(combined, motif_id, p_combined)
  combined[]
}

# ---- internal: empty result schema ---------------------------------------

.empty_assoc_dt <- function(method) {
  data.table::data.table(
    motif_id       = character(),
    orbit_id       = integer(),
    gene           = character(),
    n_in           = integer(),
    n_out          = integer(),
    n_in_expressed = integer(),
    mean_in        = numeric(),
    mean_out       = numeric(),
    log2FC         = numeric(),
    U              = numeric(),
    z              = numeric(),
    p              = numeric(),
    p_adj_BH       = numeric()
  )
}

.empty_assoc_per_sample_dt <- function() {
  data.table::data.table(
    motif_id         = character(),
    orbit_id         = integer(),
    gene             = character(),
    n_samples_tested = integer(),
    n_in_total       = integer(),
    n_out_total      = integer(),
    mean_log2FC      = numeric(),
    sd_log2FC        = numeric(),
    z_combined       = numeric(),
    p_combined       = numeric(),
    p_adj_BH         = numeric()
  )
}
