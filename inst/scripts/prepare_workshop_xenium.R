#!/usr/bin/env Rscript

# Prepare the workshop Xenium dataset for smotif.
#
# Stage A (slow, cached):  build an annotated Giotto object
#   - import via createGiottoXeniumObject
#   - QC, normalize, PCA, expression kNN, Leiden at res 0.4 and 0.8
#   - top-15 markers per cluster (gini), one-vs-all
#   - LLM cell-type calls (or fallback to cluster_N labels), with caching
#   - persist via saveGiotto
#
# Stage B (fast, recompute cheap): derive parquet + niches
#   - extract nodes (cell_id, x, y, z, cell_type, sample_id, region_id)
#   - niche clustering (k = 3) over k=10 spatial-neighbor cell-type composition
#   - analytical kNN graph (k = 6) for downstream motif finding
#   - write nodes.parquet, edges.parquet, expression.rds, README.md
#
# Idempotent. Re-run with caches present is fast (stage B only).

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(Matrix)
})

# Disable Giotto's Python/conda checks. None of the routines we call need
# Python, but Giotto initializes the path eagerly and crashes when no conda
# is installed. Tested against Giotto 4.2.3 / GiottoClass 0.5.1.
options(giotto.use_conda = FALSE,
        giotto.has_python = FALSE,
        giotto.no_python_warn = TRUE)
Sys.setenv(GIOTTO_PYTHON = "")

# ---- 0. paths & config -------------------------------------------------------

.locate_self <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) {
    return(normalizePath(sub("^--file=", "", m[1])))
  }
  # fallback when sourced interactively
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(f)) {
    return(normalizePath(f))
  }
  stop("Run this script with Rscript so its location is known.")
}

.find_pkg_root <- function(start) {
  d <- dirname(start)
  while (!file.exists(file.path(d, "DESCRIPTION")) && d != dirname(d)) {
    d <- dirname(d)
  }
  if (!file.exists(file.path(d, "DESCRIPTION"))) {
    stop("Could not locate the smotif package root (no DESCRIPTION found).")
  }
  d
}

PKG_ROOT <- .find_pkg_root(.locate_self())
EXTDATA  <- file.path(PKG_ROOT, "inst", "extdata", "workshop_xenium")
RAW_DIR  <- "/Users/rubendries/Documents/Datasets/Xenium/workshop_xenium"
SAMPLE_ID <- "workshop_xenium"

dir.create(EXTDATA, recursive = TRUE, showWarnings = FALSE)

# Source small helpers from the package source tree (the prep script runs
# *before* the user has prepared its data, so it cannot rely on the package
# being installed).
source(file.path(PKG_ROOT, "R", "utils.R"))

# Tunables
LEIDEN_RES_PRIMARY    <- 0.4
LEIDEN_RES_SECONDARY  <- 0.8
N_TOP_MARKERS         <- 15L
SPATIAL_K_NICHE       <- 10L     # neighbors used for niche composition
SPATIAL_K_ANALYTICAL  <- 6L      # neighbors for the motif graph
N_NICHES              <- 3L
QV_THRESHOLD          <- 20
LLM_MODEL             <- Sys.getenv("SMOTIF_LLM_MODEL", "claude-opus-4-7")
ANTHROPIC_VERSION     <- "2023-06-01"

ANNOTATED_DIR    <- file.path(EXTDATA, "giotto_annotated")
TISSUE_CTX_FILE  <- file.path(EXTDATA, "tissue_context.txt")
LLM_CACHE_FILE   <- file.path(EXTDATA, "llm_celltype_calls.json")
MAPPING_FILE     <- file.path(EXTDATA, "celltype_mapping.csv")
CLUSTER_MARK_FILE<- file.path(EXTDATA, "cluster_markers.csv")
NICHE_CENT_FILE  <- file.path(EXTDATA, "niche_centroids.csv")
NODES_FILE       <- file.path(EXTDATA, "nodes.parquet")
EDGES_FILE       <- file.path(EXTDATA, "edges.parquet")
EXPR_FILE        <- file.path(EXTDATA, "expression.rds")
README_FILE      <- file.path(EXTDATA, "README.md")

# ---- 1. tissue context -------------------------------------------------------
# The dataset's experiment.xenium identifies it as human lung cancer FFPE.
# Auto-write tissue_context.txt unless the user already provided one.
if (!file.exists(TISSUE_CTX_FILE)) {
  writeLines(
    "human lung cancer FFPE tissue (Xenium Human Multi-Tissue and Cancer panel, 397 genes)",
    TISSUE_CTX_FILE
  )
}
TISSUE_CONTEXT <- tryCatch(
  paste(readLines(TISSUE_CTX_FILE, warn = FALSE), collapse = " "),
  error = function(e) ""
)

# ---- 2. helpers --------------------------------------------------------------

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("'%s' is required for this script. Install it first.", pkg))
  }
}

# Attach a package onto the search path. Some Giotto routines (notably
# `update_giotto_params`) resolve the name of their caller via `get(name,
# envir = ...)` and only succeed if the package is attached, not merely
# loaded. Use this for Giotto and GiottoClass; namespace-based access via
# `pkg::fn` is enough for everything else.
attach_pkg <- function(pkg) {
  require_pkg(pkg)
  if (!paste0("package:", pkg) %in% search()) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE, quietly = TRUE)
    )
  }
}

# Build the prompt sent to Claude. Each cluster appears as a section listing
# its top markers and (when available) effect sizes.
.build_llm_prompt <- function(cluster_markers, tissue_context) {
  header <- if (nzchar(tissue_context)) {
    sprintf("Tissue context: %s.\n\n", tissue_context)
  } else {
    ""
  }
  sections <- vapply(seq_along(cluster_markers), function(i) {
    cl <- names(cluster_markers)[i]
    df <- cluster_markers[[i]]
    lines <- if (!is.null(df$score)) {
      paste0("  - ", df$gene, " (", signif(df$score, 3), ")")
    } else {
      paste0("  - ", df$gene)
    }
    paste0("Cluster ", cl, ":\n", paste(lines, collapse = "\n"))
  }, character(1))
  body <- paste(sections, collapse = "\n\n")
  paste0(
    header,
    "You are a spatial-transcriptomics annotator. For each cluster below, ",
    "infer the most likely cell type from the listed top marker genes ",
    "(ranked by gini effect size when shown in parentheses). Use the tissue ",
    "context to disambiguate when multiple cell types are plausible.\n\n",
    body,
    "\n\nReturn ONLY a JSON array, no prose, with one object per cluster:\n",
    '[{"cluster":"<id>","cell_type":"<short label>",',
    '"confidence":"high|medium|low","rationale":"<one short sentence naming key markers>"}]\n',
    "Use snake_case_lowercase cell_type labels (e.g. 'tumor_epithelial', ",
    "'cd8_t_cell', 'macrophage', 'fibroblast'). Keep rationale under 25 words."
  )
}

.call_anthropic <- function(prompt, model, key) {
  require_pkg("httr2")
  require_pkg("jsonlite")
  resp <- httr2::request("https://api.anthropic.com/v1/messages") |>
    httr2::req_headers(
      "x-api-key" = key,
      "anthropic-version" = ANTHROPIC_VERSION,
      "content-type" = "application/json"
    ) |>
    httr2::req_body_json(list(
      model = model,
      max_tokens = 4000,
      messages = list(list(role = "user", content = prompt))
    )) |>
    httr2::req_timeout(120) |>
    httr2::req_perform()
  parsed <- httr2::resp_body_json(resp)
  text <- parsed$content[[1]]$text
  # Strip optional ```json fences if the model returns them.
  text <- sub("^\\s*```(?:json)?\\s*", "", text)
  text <- sub("\\s*```\\s*$", "", text)
  jsonlite::fromJSON(text)
}

.validate_calls <- function(calls, expected_clusters) {
  stopifnot(is.data.frame(calls))
  needed <- c("cluster", "cell_type", "confidence", "rationale")
  missing_cols <- setdiff(needed, colnames(calls))
  if (length(missing_cols)) {
    stop("LLM JSON missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  calls$cluster <- as.character(calls$cluster)
  miss_cl <- setdiff(as.character(expected_clusters), calls$cluster)
  if (length(miss_cl)) {
    stop("LLM did not return calls for clusters: ", paste(miss_cl, collapse = ", "))
  }
  bad_conf <- !calls$confidence %in% c("high", "medium", "low")
  if (any(bad_conf)) {
    warning(
      sum(bad_conf), " call(s) had non-standard confidence values; coercing to 'low'."
    )
    calls$confidence[bad_conf] <- "low"
  }
  if (any(is.na(calls$cell_type) | !nzchar(calls$cell_type))) {
    stop("LLM returned an empty cell_type for at least one cluster.")
  }
  calls
}

# Compute, for each cell, the cell-type composition of its k spatial nearest
# neighbors (including self). Uses only x,y; z is ignored when missing.
.compose_neighborhoods <- function(coords, cell_types, k) {
  require_pkg("FNN")
  nn <- FNN::get.knn(coords, k = k - 1L)
  # idx is k-1 NNs; prepend self so we get k cells per neighborhood.
  N <- nrow(coords)
  idx <- cbind(seq_len(N), nn$nn.index)
  ct <- factor(cell_types)
  L <- nlevels(ct)
  # For each row, count neighbors per type, then divide by k.
  vec <- matrix(as.integer(ct)[as.vector(idx)], nrow = N, ncol = k)
  comp <- matrix(0.0, nrow = N, ncol = L,
                 dimnames = list(NULL, levels(ct)))
  for (j in seq_len(L)) {
    comp[, j] <- rowSums(vec == j) / k
  }
  comp
}

# kNN edge list (k neighbors per node), undirected, deduped, source < target.
.knn_edges <- function(coords, cell_ids, k, max_distance = NULL) {
  require_pkg("FNN")
  nn <- FNN::get.knn(coords, k = k)
  N  <- nrow(coords)
  src <- rep(seq_len(N), each = k)
  tgt <- as.vector(t(nn$nn.index))
  dst <- as.vector(t(nn$nn.dist))
  # canonicalize: source < target index, drop self-edges (none from get.knn)
  swap <- src > tgt
  s2 <- ifelse(swap, tgt, src)
  t2 <- ifelse(swap, src, tgt)
  dt <- data.table::data.table(s = s2, t = t2, distance = dst)
  if (!is.null(max_distance)) dt <- dt[distance <= max_distance]
  dt <- unique(dt, by = c("s", "t"))
  data.table::data.table(
    source    = cell_ids[dt$s],
    target    = cell_ids[dt$t],
    distance  = dt$distance,
    sample_id = SAMPLE_ID
  )
}

# ---- 3. STAGE A: annotated Giotto object ------------------------------------

run_stage_a <- function() {
  attach_pkg("Giotto")
  attach_pkg("GiottoClass")
  .step_log("stage A: importing Xenium dataset")
  gobj <- Giotto::createGiottoXeniumObject(
    xenium_dir          = RAW_DIR,
    qv_threshold        = QV_THRESHOLD,
    load_images         = NULL,
    load_aligned_images = NULL,
    load_transcripts    = FALSE,
    load_expression     = TRUE,
    load_cellmeta       = TRUE,
    verbose             = FALSE
  )

  .step_log("stage A: QC + normalize")
  # Workaround: Giotto's update_giotto_params resolves the names of caller
  # and helper functions via get(name, envir = ...) and won't find unexported
  # GiottoClass internals through the search path. Hoist `subsetGiotto` into
  # globalenv. Tested against Giotto 4.2.3 / GiottoClass 0.5.1.
  if (!exists("subsetGiotto", envir = globalenv(), inherits = FALSE)) {
    assign("subsetGiotto",
           getFromNamespace("subsetGiotto", "GiottoClass"),
           envir = globalenv())
  }
  gobj <- filterGiotto(
    gobject = gobj,
    expression_threshold = 1,
    feat_det_in_min_cells = 3,
    min_det_feats_per_cell = 5,
    verbose = FALSE
  )
  gobj <- Giotto::normalizeGiotto(gobj, verbose = FALSE)

  .step_log("stage A: PCA + expression kNN")
  # Small panel (397 genes): skip HVG selection and use all features.
  gobj <- Giotto::runPCA(gobj, feats_to_use = NULL, verbose = FALSE)
  gobj <- Giotto::createNearestNetwork(
    gobj, type = "sNN", dim_reduction_to_use = "pca",
    dimensions_to_use = 1:20, k = 30, verbose = FALSE
  )

  .step_log("stage A: Leiden clustering at two resolutions")
  gobj <- Giotto::doLeidenCluster(
    gobj, resolution = LEIDEN_RES_PRIMARY,
    name = "leiden_0.4", seed_number = 1234
  )
  gobj <- Giotto::doLeidenCluster(
    gobj, resolution = LEIDEN_RES_SECONDARY,
    name = "leiden_0.8", seed_number = 1234
  )

  .step_log("stage A: spatial kNN network for niche composition")
  gobj <- Giotto::createSpatialNetwork(
    gobj, method = "kNN", k = SPATIAL_K_NICHE,
    name = "spatial_knn_niche", return_gobject = TRUE, verbose = FALSE
  )

  .step_log("stage A: gini markers, one-vs-all on leiden_0.4")
  markers_dt <- Giotto::findGiniMarkers_one_vs_all(
    gobject = gobj,
    cluster_column = "leiden_0.4",
    expression_values = "normalized",
    verbose = FALSE
  )
  data.table::setDT(markers_dt)
  # The score column is named differently across Giotto versions; resolve it.
  cand <- intersect(c("comb_score", "expression_gini", "score"),
                    colnames(markers_dt))
  score_col <- if (length(cand)) cand[1] else NA_character_
  gene_col <- intersect(c("feats", "feature", "gene", "gene_id"),
                        colnames(markers_dt))[1]
  cluster_col <- intersect(c("cluster", "cluster_id"), colnames(markers_dt))[1]
  if (is.na(gene_col) || is.na(cluster_col)) {
    stop("Could not detect gene/cluster columns from findGiniMarkers_one_vs_all output: ",
         paste(colnames(markers_dt), collapse = ", "))
  }
  data.table::setnames(markers_dt, gene_col, "gene")
  data.table::setnames(markers_dt, cluster_col, "cluster")
  if (!is.na(score_col)) data.table::setnames(markers_dt, score_col, "score")
  markers_dt[, cluster := as.character(cluster)]

  cluster_ids <- sort(unique(markers_dt$cluster))
  per_cluster <- lapply(cluster_ids, function(cl) {
    sub <- markers_dt[cluster == cl]
    if (!is.null(sub$score)) data.table::setorder(sub, -score)
    head(sub, N_TOP_MARKERS)
  })
  names(per_cluster) <- cluster_ids
  marker_hash <- .smotif_hash(list(
    cluster_ids,
    lapply(per_cluster, function(d) sort(d$gene))
  ))

  # Resolve cell-type labels: hand-curated CSV (newest) > LLM cache (matching
  # hash) > fresh LLM call > cluster fallback.
  call_source <- "unset"
  calls <- NULL

  if (file.exists(MAPPING_FILE) && file.exists(LLM_CACHE_FILE) &&
      file.info(MAPPING_FILE)$mtime >= file.info(LLM_CACHE_FILE)$mtime) {
    .step_log("stage A: using hand-curated celltype_mapping.csv")
    calls <- data.table::fread(MAPPING_FILE, colClasses = list(character = "cluster"))
    call_source <- "human_curated_csv"
  } else if (file.exists(LLM_CACHE_FILE)) {
    cached <- jsonlite::fromJSON(LLM_CACHE_FILE, simplifyVector = TRUE)
    if (!is.null(cached$marker_hash) && identical(cached$marker_hash, marker_hash)) {
      .step_log("stage A: reusing cached LLM cell-type calls (hash match)")
      calls <- as.data.table(cached$calls)
      call_source <- sprintf("llm_cache (model=%s)", cached$model %||% "?")
    } else {
      .step_log("stage A: LLM cache exists but marker hash changed; will re-call")
    }
  }

  if (is.null(calls)) {
    key <- Sys.getenv("ANTHROPIC_API_KEY")
    if (nzchar(key)) {
      .step_log(sprintf("stage A: calling Anthropic API (%s)", LLM_MODEL))
      cluster_markers_for_prompt <- lapply(per_cluster, function(d) {
        list(gene = d$gene, score = if (!is.null(d$score)) d$score else NULL)
      })
      prompt <- .build_llm_prompt(cluster_markers_for_prompt, TISSUE_CONTEXT)
      raw <- .call_anthropic(prompt, model = LLM_MODEL, key = key)
      calls <- as.data.table(.validate_calls(raw, cluster_ids))
      jsonlite::write_json(
        list(
          marker_hash = marker_hash,
          model       = LLM_MODEL,
          tissue_context = TISSUE_CONTEXT,
          generated_at   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
          calls       = calls
        ),
        LLM_CACHE_FILE, pretty = TRUE, auto_unbox = TRUE
      )
      data.table::fwrite(calls, MAPPING_FILE)
      call_source <- sprintf("llm_fresh (model=%s)", LLM_MODEL)
    } else {
      warning(
        "ANTHROPIC_API_KEY is unset. Falling back to cluster_<n> labels. ",
        "Set the key and rerun, or hand-curate ",
        basename(MAPPING_FILE), " (cluster, cell_type, confidence, rationale) ",
        "and rerun.", call. = FALSE, immediate. = TRUE
      )
      # Persist top markers per cluster so the user has something to read.
      mark_long <- data.table::rbindlist(lapply(cluster_ids, function(cl) {
        d <- per_cluster[[cl]]
        data.table::data.table(
          cluster = cl,
          rank    = seq_len(nrow(d)),
          gene    = d$gene,
          score   = if (!is.null(d$score)) d$score else NA_real_
        )
      }))
      data.table::fwrite(mark_long, CLUSTER_MARK_FILE)
      calls <- data.table::data.table(
        cluster   = cluster_ids,
        cell_type = paste0("cluster_", cluster_ids),
        confidence = "low",
        rationale = "fallback: LLM not called (ANTHROPIC_API_KEY unset)"
      )
      call_source <- "cluster_fallback"
    }
  }

  # Apply mapping: every cell inherits its cluster's cell_type label.
  meta <- Giotto::pDataDT(gobj)
  meta_dt <- as.data.table(meta)[, .(cell_ID, leiden_0.4 = as.character(leiden_0.4))]
  meta_dt <- merge(
    meta_dt, calls[, .(cluster, cell_type)],
    by.x = "leiden_0.4", by.y = "cluster", all.x = TRUE, sort = FALSE
  )
  data.table::setkey(meta_dt, cell_ID)
  gobj <- Giotto::addCellMetadata(
    gobj,
    new_metadata = meta_dt[, .(cell_ID, cell_type)],
    by_column = TRUE,
    column_cell_ID = "cell_ID"
  )

  .step_log("stage A: saveGiotto")
  if (dir.exists(ANNOTATED_DIR)) unlink(ANNOTATED_DIR, recursive = TRUE)
  Giotto::saveGiotto(
    gobj,
    foldername = "giotto_annotated",
    dir = EXTDATA,
    method = "RDS",
    overwrite = TRUE,
    verbose = FALSE
  )

  list(gobj = gobj, calls = calls, call_source = call_source)
}

# ---- 4. STAGE B: derived parquet + niches -----------------------------------

run_stage_b <- function(stage_a) {
  attach_pkg("Giotto")
  attach_pkg("GiottoClass")
  gobj <- stage_a$gobj
  .step_log("stage B: extracting cells and coordinates")

  meta <- as.data.table(Giotto::pDataDT(gobj))
  spat <- as.data.table(Giotto::getSpatialLocations(gobj, output = "data.table"))
  # Spatial location columns are typically sdimx, sdimy (and sdimz if 3D).
  xcol <- intersect(c("sdimx", "x", "x_centroid"), colnames(spat))[1]
  ycol <- intersect(c("sdimy", "y", "y_centroid"), colnames(spat))[1]
  zcol <- intersect(c("sdimz", "z", "z_centroid"), colnames(spat))
  cell_id_col <- intersect(c("cell_ID", "cell_id"), colnames(spat))[1]
  if (is.na(xcol) || is.na(ycol) || is.na(cell_id_col)) {
    stop("Could not detect coords or cell-id columns in spatial locations: ",
         paste(colnames(spat), collapse = ", "))
  }
  spat <- spat[, c(cell_id_col, xcol, ycol, zcol), with = FALSE]
  data.table::setnames(spat, c(cell_id_col, xcol, ycol),
                       c("cell_id", "x", "y"))
  if (length(zcol) == 1L) {
    data.table::setnames(spat, zcol, "z")
  } else {
    spat[, z := NA_real_]
  }
  data.table::setkey(spat, cell_id)

  meta <- meta[, .(cell_id = cell_ID,
                   cell_type = as.character(cell_type),
                   leiden_0.4 = as.character(leiden_0.4),
                   leiden_0.8 = as.character(leiden_0.8))]
  nodes <- merge(spat, meta, by = "cell_id", sort = FALSE)
  nodes[, sample_id := SAMPLE_ID]

  if (anyNA(nodes$cell_type)) {
    n_na <- sum(is.na(nodes$cell_type))
    warning(n_na, " cell(s) lack a cell_type label after stage A; coercing to 'unknown'.")
    nodes[is.na(cell_type), cell_type := "unknown"]
  }

  .step_log("stage B: niche clustering (k = 3) over cell-type composition")
  coords <- as.matrix(nodes[, .(x, y)])
  comp <- .compose_neighborhoods(coords, nodes$cell_type, k = SPATIAL_K_NICHE)
  set.seed(1)
  km <- kmeans(comp, centers = N_NICHES, nstart = 25, iter.max = 50)
  nodes[, region_id := paste0("niche_", km$cluster)]

  niche_centroids <- as.data.table(km$centers, keep.rownames = "niche")
  niche_centroids[, niche := paste0("niche_", niche)]
  data.table::fwrite(niche_centroids, NICHE_CENT_FILE)

  .step_log("stage B: building analytical kNN graph (k = 6)")
  edges <- .knn_edges(coords, nodes$cell_id, k = SPATIAL_K_ANALYTICAL)

  # Final node columns in the order specified by the data model.
  nodes_out <- nodes[, .(cell_id, x, y, z, cell_type, sample_id, region_id,
                         leiden_0.4, leiden_0.8)]

  .step_log("stage B: writing nodes.parquet, edges.parquet, expression.rds")
  arrow::write_parquet(nodes_out, NODES_FILE)
  arrow::write_parquet(edges,     EDGES_FILE)

  expr <- Giotto::getExpression(
    gobj, values = "raw", spat_unit = "cell", feat_type = "rna",
    output = "matrix"
  )
  if (!inherits(expr, "dgCMatrix")) expr <- as(expr, "CsparseMatrix")
  saveRDS(expr, EXPR_FILE, compress = "xz")

  list(
    nodes = nodes_out,
    edges = edges,
    niche_centroids = niche_centroids,
    expr = expr,
    call_source = stage_a$call_source,
    calls = stage_a$calls
  )
}

# ---- 5. README ---------------------------------------------------------------

write_readme <- function(stage_b) {
  nodes <- stage_b$nodes
  edges <- stage_b$edges
  cn    <- stage_b$niche_centroids
  calls <- stage_b$calls

  ct_tab <- nodes[, .N, by = cell_type][order(-N)]
  rg_tab <- nodes[, .N, by = region_id][order(region_id)]
  ed_dist <- if (nrow(edges)) {
    quantile(edges$distance, c(0, 0.5, 0.95, 1.0), na.rm = TRUE)
  } else {
    c(NA, NA, NA, NA)
  }

  size4_note <- if (nrow(nodes) > 50000L) {
    paste0("- Cell count exceeds 50k; size-4 motif enumeration should use ",
           "`anchored_on` sampling.\n")
  } else {
    paste0("- Cell count is ", nrow(nodes), " (< 50k); size-4 enumeration ",
           "should run on the full graph without `anchored_on` sampling.\n")
  }

  ctx_note <- if (file.exists(TISSUE_CTX_FILE)) {
    sprintf("Tissue context: %s\n", TISSUE_CONTEXT)
  } else {
    "Tissue context: NOT SET (LLM was prompted without context).\n"
  }

  fallback_note <- switch(
    sub(" .*", "", stage_b$call_source),
    "cluster_fallback"   = paste0("**LLM fallback fired:** `ANTHROPIC_API_KEY` was unset, ",
                                  "so cells are labeled `cluster_<n>`. To get real cell-type ",
                                  "names, either (a) set the API key and rerun, or (b) ",
                                  "hand-edit `celltype_mapping.csv` (see `cluster_markers.csv` ",
                                  "for top markers per cluster) and rerun.\n"),
    "human_curated_csv"  = "Cell-type labels were taken from a hand-curated `celltype_mapping.csv`.\n",
    "llm_cache"          = paste0("Cell-type labels reused from cached LLM call ",
                                  "(`llm_celltype_calls.json`).\n"),
    "llm_fresh"          = "Cell-type labels obtained from a fresh LLM call.\n",
    sprintf("Cell-type call source: %s.\n", stage_b$call_source)
  )

  low_conf <- if (!is.null(calls$confidence)) {
    sum(tolower(calls$confidence) == "low")
  } else 0L
  conf_note <- if (low_conf > 0L) {
    sprintf("- %d cluster(s) flagged as low confidence — review `celltype_mapping.csv`.\n",
            low_conf)
  } else ""

  fmt_table <- function(dt) {
    if (!nrow(dt)) return("")
    cols <- colnames(dt)
    head <- paste0("| ", paste(cols, collapse = " | "), " |")
    sep  <- paste0("|", paste(rep("---", length(cols)), collapse = "|"), "|")
    rows <- vapply(seq_len(nrow(dt)), function(i) {
      paste0("| ", paste(as.character(unlist(dt[i])), collapse = " | "), " |")
    }, character(1))
    paste(c(head, sep, rows), collapse = "\n")
  }

  rationale_table <- if (!is.null(calls$rationale)) {
    fmt_table(calls[, .(cluster, cell_type, confidence, rationale)])
  } else ""

  body <- paste0(
    "# workshop_xenium — prepared for `smotif`\n\n",
    "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n",
    "Source: `", RAW_DIR, "`\n\n",
    "## Stage A summary\n\n",
    ctx_note,
    fallback_note,
    "\n",
    "Cell-type call rationales:\n\n",
    rationale_table, "\n\n",
    conf_note,
    "## Cells\n\n",
    "- Total cells (after QC): **", nrow(nodes), "**\n",
    size4_note,
    "\nCell-type distribution:\n\n",
    fmt_table(ct_tab), "\n\n",
    "## Niches\n\n",
    "Niche assignments (k-means, k = ", N_NICHES, ", on cell-type composition over ",
    SPATIAL_K_NICHE, " spatial neighbors + self):\n\n",
    fmt_table(rg_tab), "\n\n",
    "Niche centroid composition (rows sum to 1):\n\n",
    fmt_table(cn), "\n\n",
    "## Edges (analytical kNN, k = ", SPATIAL_K_ANALYTICAL, ", undirected)\n\n",
    "- Edge count: **", nrow(edges), "**\n",
    "- Distance: min ", signif(ed_dist[1], 4), ", median ", signif(ed_dist[2], 4),
    ", p95 ", signif(ed_dist[3], 4), ", max ", signif(ed_dist[4], 4), "\n\n",
    "## What to verify\n\n",
    "- [ ] Cell-type labels look biologically plausible for human lung cancer FFPE.\n",
    "- [ ] Each niche has a non-trivial number of cells (no single tiny niche).\n",
    "- [ ] Niche centroid compositions are interpretable (e.g. one tumor-dense, ",
    "one stromal/immune, one mixed) — this is what motif tests will be ",
    "stratified on.\n",
    "- [ ] Edge-length p95 looks reasonable for the tissue scale (a long ",
    "p95 hints that some cells are isolated and pulling in distant kNN).\n"
  )
  writeLines(body, README_FILE)
  invisible(README_FILE)
}

# ---- 6. main -----------------------------------------------------------------

main <- function() {
  if (dir.exists(ANNOTATED_DIR)) {
    attach_pkg("Giotto"); attach_pkg("GiottoClass")
    .step_log("stage A: loading cached annotated Giotto object")
    # python_path = NA disables the Python/conda check we don't need; this
    # script never touches the Python-backed routines.
    gobj <- GiottoClass::loadGiotto(
      path_to_folder = ANNOTATED_DIR,
      reconnect_giottoImage = FALSE,
      init_gobject = FALSE,
      python_path = NA,
      verbose = FALSE
    )
    # Re-derive call_source from on-disk files for the README.
    call_source <- if (file.exists(LLM_CACHE_FILE)) "llm_cache (cached)" else "cluster_fallback (cached)"
    calls <- if (file.exists(MAPPING_FILE)) {
      data.table::fread(MAPPING_FILE, colClasses = list(character = "cluster"))
    } else NULL
    stage_a <- list(gobj = gobj, calls = calls, call_source = call_source)
  } else {
    stage_a <- run_stage_a()
  }
  stage_b <- run_stage_b(stage_a)
  write_readme(stage_b)
  .step_log(sprintf("done. README written to %s", README_FILE))
  invisible(NULL)
}

main()
