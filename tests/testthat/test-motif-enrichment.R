# `motif_enrichment()` has two backends that must not drift apart. The Rust
# path derives canonical colored forms from each topology's automorphism group;
# the pure-R path hard-codes the same forms for the two size-3 topologies. The
# cross-backend test is what keeps those two descriptions of the same thing in
# agreement -- it is the reason the R path is worth having at all.

.sg <- function(n = 200, k = 5, types = c("T", "B", "M"), seed = 3) {
    set.seed(seed)
    build_spatial_graph(
        cbind(x = runif(n), y = runif(n)),
        cell_types = sample(types, n, TRUE),
        sample_id = "s1", method = "knn", k = k
    )
}

test_that("the two backends agree exactly at sizes 2 and 3", {
    skip_if_not_installed("smotifrs")
    sg <- .sg()
    for (sz in 2:3) {
        a <- motif_enrichment(sg, size = sz, n_perm = 100L, seed = 9L, backend = "rust")
        b <- motif_enrichment(sg, size = sz, n_perm = 100L, seed = 9L, backend = "r")
        expect_identical(
            attr(a, "n_instances"), attr(b, "n_instances"),
            info = paste("size", sz)
        )
        expect_setequal(a$motif_id, b$motif_id)
        m <- merge(
            a[, list(motif_id, ru = observed)],
            b[, list(motif_id, rr = observed)],
            by = "motif_id"
        )
        expect_equal(m$ru, m$rr, info = paste("size", sz))
    }
})

test_that("observed counts sum to the instance count in both backends", {
    sg <- .sg(n = 120)
    b <- motif_enrichment(sg, size = 3L, n_perm = 50L, seed = 1L, backend = "r")
    expect_equal(sum(b$observed), attr(b, "n_instances"))
    skip_if_not_installed("smotifrs")
    a <- motif_enrichment(sg, size = 3L, n_perm = 50L, seed = 1L, backend = "rust")
    expect_equal(sum(a$observed), attr(a, "n_instances"))
})

test_that("the result contract is stable across backends", {
    sg <- .sg(n = 120)
    want <- c(
        "motif_id", "topology", "size", "color_tuple", "observed", "expected",
        "sd_null", "z", "fold", "p_enrich", "p_deplete", "p_adj"
    )
    b <- motif_enrichment(sg, size = 3L, n_perm = 50L, seed = 1L, backend = "r")
    expect_identical(names(b), want)
    expect_type(b$color_tuple, "list")
    expect_true(all(b$p_enrich > 0 & b$p_enrich <= 1))
    skip_if_not_installed("smotifrs")
    a <- motif_enrichment(sg, size = 3L, n_perm = 50L, seed = 1L, backend = "rust")
    expect_identical(names(a), want)
})

test_that("size 4 in pure R is refused with a reason, not attempted", {
    sg <- .sg(n = 80)
    expect_error(
        motif_enrichment(sg, size = 4L, backend = "r"),
        "needs the smotifrs backend"
    )
})

test_that("the conditional null is refused by the pure-R backend", {
    sg <- .sg(n = 80)
    expect_error(
        motif_enrichment(sg, size = 3L, backend = "r", null = "conditional"),
        "smotifrs backend only"
    )
})

test_that("anchoring restricts to subgraphs touching the anchor set", {
    sg <- .sg(n = 120)
    anchors <- sg$nodes$cell_id[1:20]
    full <- motif_enrichment(sg, size = 3L, n_perm = 20L, seed = 1L, backend = "r")
    anc <- motif_enrichment(sg,
        size = 3L, n_perm = 20L, seed = 1L,
        backend = "r", anchored_on = anchors
    )
    expect_lt(attr(anc, "n_instances"), attr(full, "n_instances"))
    expect_gt(attr(anc, "n_instances"), 0)
})

test_that("a bad anchor is reported rather than silently dropped", {
    sg <- .sg(n = 60)
    expect_error(
        motif_enrichment(sg, size = 3L, backend = "r", anchored_on = "not_a_cell"),
        "not present"
    )
})

test_that("stratifying needs a column that exists", {
    sg <- .sg(n = 60)
    expect_error(
        motif_enrichment(sg, size = 3L, backend = "r", null = "stratified"),
        "region_id"
    )
    expect_error(
        motif_enrichment(sg, size = 3L, backend = "r", strata_column = "nope"),
        "not found"
    )
})

test_that("results are reproducible for a given seed", {
    sg <- .sg(n = 100)
    a <- motif_enrichment(sg, size = 3L, n_perm = 40L, seed = 5L, backend = "r")
    b <- motif_enrichment(sg, size = 3L, n_perm = 40L, seed = 5L, backend = "r")
    expect_equal(a$expected, b$expected)
})
