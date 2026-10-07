test_that("motif_enrichment_stream delegates to the backend", {
  skip_if_not_installed("smotifrs")
  skip_if_not("motif_enrichment_stream" %in% getNamespaceExports("smotifrs"))
  skip_if_not_installed("arrow")
  skip_if_not_installed("nanoarrow")
  set.seed(1)
  n <- 60L
  el <- igraph::as_edgelist(igraph::sample_gnp(n, 0.1), names = FALSE)
  ct <- factor(sample(c("A", "B"), n, TRUE))
  edges <- function() {
    arrow::arrow_table(from_id = as.integer(el[, 1]), to_id = as.integer(el[, 2]))
  }
  a <- motif_enrichment_stream(edges(), seq_len(n), ct,
                               size = 3L, n_perm = 19L, seed = 2L)
  b <- smotifrs::motif_enrichment_stream(edges(), seq_len(n), ct,
                                         size = 3L, n_perm = 19L, seed = 2L)
  expect_identical(a, b)
})
