
# RIPPLET

<!-- badges: start -->

<!-- add CI/coverage/pkgdown badges here when available -->

<!-- badges: end -->

RIPPLET (**R**andom-walk-based **I**nference of Gene & **P**athway
**P**erturbations integrating Sample-**L**evel and Cohort-wide
**E**ffec**T**s) is a DNA-only framework for sample-level driver and
pathway perturbation scoring. It:

- builds and normalizes a network adjacency matrix from an edgelist,
- propagates binary mutations via Random Walk with Restart (RWR),
- derives sample–sample similarity in a low-dimensional SVD space, and
- computes gene impact scores by weighting each sample-gene score
  according to their cohort neighbors.

Optionally, RIPPLET projects GIS onto topological pathway graphs
(SPIA-style) to yield normalized pathway perturbation scores (NPPS).

RIPPLET is designed for robust, transcriptome-independent driver gene &
pathway inference in large cohorts and clinical sequencing contexts.

## Installation

You can install the development version of RIPPLET from GitHub:

``` r
install.packages("devtools")
devtools::install_github("wiesmannf/RIPPLET")
```

## Quick start: Gene Impact Scores (recommended)

Returns gene impact scores to be used in driver prediction and for
pathway perturbation calculation.

``` r
library(RIPPLET)

# Demo inputs shipped with the package (derived from TCGA-SKCM)
data(edgelist) # example STRING-like network (tissue-contextualized for SKCM)
data(mutmat)   # example TCGA-SKCM 0/1 mutation matrix (genes x samples)

# End-to-end gene pipeline: network → normalize → RWR → similarity → GIS
gis <- ripplet_gene(
  edgelist = edgelist,
  mutmat   = mutmat,
  alpha    = 0.7,     # RWR restart probability
  maxiters = 20000,   # RWR max iterations
  eps      = 1e-14,   # RWR convergence tolerance
  k        = 20,      # SVD rank for sample similarity
  delta    = 0.2,     # cohort-weight decay
  verbose  = TRUE
)

# Top genes for the first sample
head(sort(gis[, 1], decreasing = TRUE), 10)
```

**What you get:** `gis` is a genes × samples matrix with cohort-aware
gene-impact scores scaled to \[0,1\] per sample.

## Gene Impact Scores: Minimal performance check (optional)

RIPPLET reports **normalized partial AUCs (npAUCs)** for precision,
recall, and F1 up to a dynamically chosen rank, following the evaluation
strategy introduced in *PersonaDrive* (Erten *et al.*, 2022,
*Bioinformatics*, doi:
[10.1093/bioinformatics/btac329](https://doi.org/10.1093/bioinformatics/btac329)).

``` r
library(RIPPLET)
data(refgenes)  # example reference gene set (SKCM)

# Suppose 'gis' is from the gene pipeline above:
perf <- personalized_perf(
  global_ref     = refgenes,
  mut_mat        = mutmat,
  impact_scores  = gis,     # gene-impact scores
  max_rank       = NULL     # adaptive rank by default
)

head(perf$avg_df)
perf$auc_stats
```

## Quick start: Pathway Perturbation Scores

**Heads-up:** Fetching & normalizing pathway graphs and running pathway
perturbation scoring can be **compute-intensive**, especially with large
pathway databases and cohorts. Consider reusing saved pathway objects
(`saveRDS(pw,"pathways.rds")`) to avoid repeated fetching.

``` r
library(RIPPLET)
data(edgelist)
data(mutmat)

# 1) Gene scores (as above)
gis <- ripplet_gene(edgelist, mutmat, verbose = FALSE)

# 2) Fetch + normalize pathway graphs and compute MaxHit denominators
#    Choose any graphite DB: "kegg", "panther", "pathbank", "pharmgkb", "reactome", "smpdb", "wikipathways"
#    See graphite::pathwayDatabases()
pw <- get_ripplet_pathways(data.base = "reactome", 
                           dampening_alpha = 0.99, 
                           show_progress = TRUE)

# 3) One-call pathway wrapper: GIS → NPPS
npps <- ripplet_path(
  gene_scores        = gis,
  pathways           = pw,        # full list from get_ripplet_pathways()
  threshold_level_ps = 0.01,      # zero-out tiny gene scores
  show_progress      = TRUE
)

dim(npps)         # pathways x samples
npps[1:5, 1:3]    # peek
```

------------------------------------------------------------------------

## Alternative — individual steps

If you prefer full control, you can run each stage explicitly.

### 1) Gene pipeline (explicit)

``` r
library(RIPPLET)
data(edgelist)
data(mutmat)

adj      <- build_network_adj(edgelist)
adj_norm <- normalize_network_adj(adj)

rwr <- rwr_multi(adj_norm, mutmat,
                 alpha    = 0.7,
                 maxiters = 20000,
                 eps      = 1e-14,
                 verbose  = FALSE)

sim <- compute_sample_similarity(rwr, k = 20, seed = 1)

gis <- impact_scoring(rwr, sim, delta = 0.2)

head(sort(gis[, 1], decreasing = TRUE), 10)
```

### 2) Pathway pipeline (explicit)

``` r
# Fetch & convert pathway graphs
pw_graphs <- get_graphs(data.base = "reactome")

# Normalize to SPIA-style adjacencies
pw_norm   <- norm_graphs(pathway_input = pw_graphs, dampening_alpha = 0.99)

# Compute per-pathway maximum propagation scores (denominator for normalization)
maxhit    <- get_pps_max_hit(pw_norm$SPIA_Norm_Adj_Matrices, show_progress = FALSE)

# Assemble full pathway object (or just keep pieces separately)
pw <- within(pw_norm, { MaxHit <- maxhit })

# Sample-specific pathway tables (PPS)
pps <- get_pps(
  gene_scores        = gis,
  pathways           = pw,        # can also pass pw$SPIA_Norm_Adj_Matrices directly
  threshold_level_ps = 0.01,
  show_progress      = TRUE
)

# Normalized pathway matrix (NPPS)
npps <- get_npps(
  pps            = pps,
  pathways       = pw,      # or use a named numeric vector of MaxHit
  show_progress  = TRUE
)
```

## Inputs & outputs at a glance

- **Network:** two-column edge list → sparse, undirected adjacency;
  row-normalized to be stochastic after `normalize_network_adj()`.
- **Mutations:** 0/1 matrix (**genes × samples**), pre-processed to
  protein-altering SNVs/InDels/CNAs/fusions per your pipeline.
- **Gene output (GIS):** numeric matrix (**genes × samples**), min–max
  scaled to **\[0,1\]** per sample.
- **Pathways (topology-aware):**
  - Pull graphs from **graphite** with `get_graphs()`.
  - Normalize to SPIA-style adjacencies with `norm_graphs()`.
  - Compute per-pathway **MaxHit** with `get_pps_max_hit()`.
  - Create per-sample pathway tables (**PPS**) with `get_pps()`.
  - Convert to normalized pathway scores (**NPPS**) by dividing by
    MaxHit via `get_npps()`.
  - Or run the one-call wrapper `ripplet_path()` (GIS → NPPS).

## Practical notes

- **Reproducibility:** `compute_sample_similarity()` seeds `irlba`
  locally (default `seed = 1`).
- **Performance:** RWR is iterative; prefer sparse matrices; a good
  starting point is `alpha ≈ 0.7`.
- **Stability:** Near-singular pathway matrices are dampened in
  `norm_graphs()` via `dampening_alpha` (default `0.99`).
- **Progress bars:** Install **progress** for nicer bars; otherwise
  messages are minimal.

## Citing RIPPLET

If you use RIPPLET, please cite the associated manuscript (details to
follow). For now, you may reference it as:

> Wiesmann F. **RIPPLET: A Mutation-Only Tool for Personalized Cancer
> Gene and Pathway Analysis in Precision Oncology.** R package version
> 0.1.0, in preparation (2025).

In addition, please also cite the following resources when applicable:

- Sales G, Calura E, Cavalieri D, Romualdi C. **graphite: a Bioconductor
  package to convert pathway topology to gene network.** *BMC
  Bioinformatics* 13, 20 (2012). doi:
  <https://doi.org/10.1186/1471-2105-13-20>  
- Han J, He Y, Li X. **PMAPscore: Identify Prognosis-Related Pathways
  Altered by Somatic Mutation.** R package version 0.1.1, 2025. CRAN:
  <https://CRAN.R-project.org/package=PMAPscore>  
- Erten M, Koyutürk M. **PersonaDrive: a personalized driver gene
  identification framework through multiomics data integration.**
  *Bioinformatics* 38(17):4095–4103 (2022). doi:
  <https://doi.org/10.1093/bioinformatics/btac329>

## Contributing / issues

Bugs, feature requests, and discussions are welcome via GitHub issues
and PRs. Please include a minimal reproducible example (data dimensions,
parameters, and session info).

## Attribution & acknowledgements

RIPPLET interfaces with graphite (AGPL-3) for pathway retrieval and ID
conversion, and adapts components from PMAPscore (GPL-2 \| GPL-3) for
pathway-level scoring. We retain upstream notices (see inst/COPYRIGHTS).
Please also cite graphite and PMAPscore.

## License

AGPL (≥ 3). See the package DESCRIPTION for details.
