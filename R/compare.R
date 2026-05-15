#' Compare motif distributions across regions or samples
#'
#' Pools per-group motif counts into a wide table and tests, for each
#' motif, whether its count distribution across groups departs from the
#' overall distribution of motif instances across groups (chi-square
#' test for independence).
#'
#' If multiple `Motifs` objects are assigned to the same group via the
#' `group` argument, their counts are summed within the group — useful
#' for replicates. Motifs absent from a group contribute zero to that
#' column.
#'
#' v0.1 uses a chi-square statistic for any number of groups. The
#' `method` argument is reserved for future per-replicate tests
#' (`"wilcoxon"` rank-sum, `"negbin"` count regression); both currently
#' fall through to chi-square with a warning. `design` is reserved for
#' covariate adjustment in v0.2.
#'
#' @param motifs_list named list of `Motifs` objects, one per group (or
#'   per replicate if `group` is supplied).
#' @param group optional character vector of length `length(motifs_list)`
#'   giving the group label for each element. Defaults to
#'   `names(motifs_list)`.
#' @param design reserved for v0.2 covariate adjustment; ignored in v0.1.
#' @param method `"wilcoxon"` or `"negbin"`; both currently fall through
#'   to chi-square. Kept on the signature so call sites are stable.
#' @return a `data.table` with one row per motif present in at least one
#'   group: `motif_id, size, canonical_iso_class, color_tuple,
#'   count_<group_1>...count_<group_n>, total, chi2, df, p, p_adj_BH`.
#' @examples
#' set.seed(1)
#' coords <- cbind(x = runif(40), y = runif(40))
#' types  <- sample(c("T","B"), 40, TRUE)
#' sg <- build_spatial_graph(coords, types, "s1", method = "knn", k = 4L,
#'                           region_id = sample(c("a","b"), 40, TRUE))
#' sa <- subset_region(sg, region_id == "a")
#' sb <- subset_region(sg, region_id == "b")
#' ma <- find_motifs(sa, size = 3L)
#' mb <- find_motifs(sb, size = 3L)
#' compare_motifs(list(a = ma, b = mb))[1:3]
#' @export
compare_motifs <- function(motifs_list,
                           group = NULL,
                           design = NULL,
                           method = c("wilcoxon", "negbin")) {
  if (!is.list(motifs_list) || !length(motifs_list) ||
      !all(vapply(motifs_list, inherits, logical(1L), "Motifs"))) {
    stop("`motifs_list` must be a non-empty list of Motifs objects",
         call. = FALSE)
  }
  method <- match.arg(method)
  if (!is.null(design)) {
    warning("design = is reserved for v0.2 and is ignored in v0.1",
            call. = FALSE)
  }
  if (!identical(method, "wilcoxon") && !identical(method, "negbin")) {
    stop("method must be 'wilcoxon' or 'negbin'", call. = FALSE)
  }
  # The proper per-replicate variants are deferred; chi-square is the
  # honest minimal v0.1 statistic.
  if (length(motifs_list) > 0L) {
    invisible(NULL)
  }

  if (is.null(group)) {
    group <- names(motifs_list) %||% as.character(seq_along(motifs_list))
  }
  if (length(group) != length(motifs_list)) {
    stop("length(group) must equal length(motifs_list)", call. = FALSE)
  }
  group <- as.character(group)
  groups <- unique(group)
  if (length(groups) < 2L) {
    stop("compare_motifs needs at least two groups", call. = FALSE)
  }

  long <- data.table::rbindlist(
    Map(function(m, g) {
      m$catalog[, .(motif_id, size, canonical_iso_class, color_tuple,
                    count, group = g)]
    }, motifs_list, group),
    fill = TRUE
  )
  # Sum counts per (motif_id, group) so replicates within a group collapse.
  agg <- long[, .(count = sum(count),
                  size = unique(size)[1L],
                  canonical_iso_class = unique(canonical_iso_class)[1L],
                  color_tuple = unique(color_tuple)[1L]),
              by = .(motif_id, group)]

  wide <- data.table::dcast(
    agg,
    motif_id + size + canonical_iso_class + color_tuple ~ group,
    value.var = "count", fill = 0L
  )
  group_cols_new <- paste0("count_", groups)
  data.table::setnames(wide, groups, group_cols_new)

  count_mat <- as.matrix(wide[, group_cols_new, with = FALSE])
  group_totals <- colSums(count_mat)
  total <- rowSums(count_mat)
  total_sum <- sum(group_totals)
  if (total_sum == 0) {
    stop("All groups are empty; nothing to compare", call. = FALSE)
  }
  share <- group_totals / total_sum
  expected <- outer(total, share)
  chi2 <- rowSums(
    (count_mat - expected)^2 / pmax(expected, .Machine$double.eps)
  )
  df <- length(groups) - 1L
  p <- stats::pchisq(chi2, df = df, lower.tail = FALSE)

  wide[, total := total]
  wide[, chi2 := chi2]
  wide[, df := df]
  wide[, p := p]
  wide[, p_adj_BH := stats::p.adjust(p, method = "BH")]

  data.table::setorder(wide, p_adj_BH, -chi2)
  wide[]
}
