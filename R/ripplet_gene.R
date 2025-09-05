#' Run the Full RIPPLET Pipeline
#'
#' From an edge list and mutation matrix, this function:
#' 1. Builds and simplifies the network adjacency matrix
#' 2. Normalizes it to row‑stochastic form
#' 3. Runs Random Walk with Restart for each sample
#' 4. Computes sample–sample similarity via truncated SVD
#' 5. Computes cohort‑informed, weighted cosine gene scores
#'
#' @param edgelist A two‑column \code{data.frame} or matrix of gene–gene edges.
#' @param mutmat A numeric matrix (genes × samples) of 0/1 mutations (or weights).
#' @param alpha Restart probability for RWR (default 0.7).
#' @param maxiters Maximum iterations for RWR (default 20000).
#' @param eps Convergence tolerance for RWR (default 1e-14).
#' @param k Number of components for sample‐similarity SVD (default 20).
#' @param delta Decay parameter for weighted scores (default 0.2).
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default \code{TRUE}).
#' @return A numeric matrix (genes × samples) of weighted cosine scores.
#' @examples
#' \dontrun{
#' # 1) build & normalize adjacency
#' adj      <- build_network_adj(edgelist)
#' adj_norm <- normalize_network_adj(adj)
#'
#' # 2) run full pipeline
#' result <- ripplet_gene(
#'   edgelist, mutmat,
#'   alpha    = 0.85,
#'   maxiters = 5000,
#'   eps      = 1e-8,
#'   k        = 10,
#'   delta    = 0.1,
#'   verbose  = FALSE
#' )
#' }
#' @export
ripplet_gene <- function(edgelist, mutmat,
                        alpha    = 0.7,
                        maxiters = 20000,
                        eps      = 1e-14,
                        k        = 20,
                        delta    = 0.2,
                        verbose  = TRUE) {
  if (verbose) message("RIPPLET [1/5]: Building network adjacency")
  network_adj      <- build_network_adj(edgelist)

  if (verbose) message("RIPPLET [2/5]: Normalizing network")
  network_adj_norm <- normalize_network_adj(network_adj)

  if (verbose) message("RIPPLET [3/5]: Running network propagation")
  rwr_profiles     <- rwr_multi(
    network_adj_norm, mutmat,
    alpha    = alpha,
    maxiters = maxiters,
    eps      = eps,
    verbose  = verbose
  )

  if (verbose) message("RIPPLET [4/5]: Computing sample similarity")
  sim_mat <- compute_sample_similarity(
    rwr_profiles,
    k    = k,
    seed = 1
  )

  if (verbose) message("RIPPLET [5/5]: Computing gene impact scores")
  result_matrix <- impact_scoring(
    rwr_profiles,
    sim_mat,
    delta = delta
  )

  if (verbose) message("RIPPLET finished")
  result_matrix
}
