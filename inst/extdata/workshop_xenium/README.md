# workshop_xenium — prepared for `smotif`

Generated: 2026-05-06 22:08:12 EDT

Source: `/Users/rubendries/Documents/Datasets/Xenium/workshop_xenium`

## Stage A summary

Tissue context: human lung cancer FFPE tissue (Xenium Human Multi-Tissue and Cancer panel, 397 genes)
**LLM fallback fired:** `ANTHROPIC_API_KEY` was unset, so cells are labeled `cluster_<n>`. To get real cell-type names, either (a) set the API key and rerun, or (b) hand-edit `celltype_mapping.csv` (see `cluster_markers.csv` for top markers per cluster) and rerun.

Cell-type call rationales:



## Cells

- Total cells (after QC): **7586**
- Cell count is 7586 (< 50k); size-4 enumeration should run on the full graph without `anchored_on` sampling.

Cell-type distribution:

| cell_type | N |
|---|---|
| cluster_4 | 1610 |
| cluster_2 | 1390 |
| cluster_1 | 1045 |
| cluster_3 | 957 |
| cluster_6 | 753 |
| cluster_5 | 725 |
| cluster_8 | 563 |
| cluster_9 | 283 |
| cluster_7 | 260 |

## Niches

Niche assignments (k-means, k = 3, on cell-type composition over 10 spatial neighbors + self):

| region_id | N |
|---|---|
| niche_1 | 4750 |
| niche_2 | 1130 |
| niche_3 | 1706 |

Niche centroid composition (rows sum to 1):

| niche | cluster_1 | cluster_2 | cluster_3 | cluster_4 | cluster_5 | cluster_6 | cluster_7 | cluster_8 | cluster_9 |
|---|---|---|---|---|---|---|---|---|---|
| niche_1 | 0.0144842105263158 | 0.284273684210529 | 0.183347368421055 | 0.0210947368421053 | 0.11642105263158 | 0.16242105263158 | 0.0560631578947364 | 0.103726315789475 | 0.0581684210526317 |
| niche_2 | 0.859203539823007 | 0.028141592920354 | 0.0336283185840709 | 0.0420353982300886 | 0.030353982300885 | 0.000353982300884956 | 0.000176991150442478 | 0.00451327433628319 | 0.0015929203539823 |
| niche_3 | 0.0174091441969519 | 0.0260257913247363 | 0.0291910902696367 | 0.840855803048066 | 0.0747362250879247 | 0.0011137162954279 | 0.00134818288393904 | 0.00861664712778428 | 0.000703399765533411 |

## Edges (analytical kNN, k = 6, undirected)

- Edge count: **26304**
- Distance: min 1.991, median 9.637, p95 19.04, max 71.05

## What to verify

- [ ] Cell-type labels look biologically plausible for human lung cancer FFPE.
- [ ] Each niche has a non-trivial number of cells (no single tiny niche).
- [ ] Niche centroid compositions are interpretable (e.g. one tumor-dense, one stromal/immune, one mixed) — this is what motif tests will be stratified on.
- [ ] Edge-length p95 looks reasonable for the tissue scale (a long p95 hints that some cells are isolated and pulling in distant kNN).

