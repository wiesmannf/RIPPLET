#' Random Walk with Restart for Multiple Samples
#'
#' Given a row‑normalized adjacency matrix and a binary mutation matrix,
#' performs a Random Walk with Restart for each patient (column of \code{mutmat})
#' until convergence or until \code{maxiters} is reached.
#'
#' @param network.adj_norm A sparse, row-normalized adjacency matrix
#'   (class \linkS4class{dgCMatrix}), e.g. from \code{\link{normalize_network_adj}}.
#' @param mutmat A numeric matrix or data‑frame with genes as rows and patients
#'   as columns; entries are 0/1. Rows (genes) not present in
#'   \code{rownames(network.adj_norm)} will be treated as zeros.
#' @param alpha Restart probability in \[0,1\]; larger means stronger pull toward
#'   initial seeds (default 0.7).
#' @param maxiters Maximum number of iterations per patient (default 20000).
#' @param eps Convergence tolerance of the score vector
#'   (default 1e-14).
#' @param verbose Logical; if \code{TRUE}, show a progress bar via
#'   \code{\link[progress]{progress_bar}} (default \code{TRUE}).
#' @return A numeric matrix (genes × patients) of steady‑state RWR scores
#'   Rows and columns carry the same dimnames as
#'   \code{network.adj_norm} and \code{mutmat}, respectively.
#' @examples
#' \dontrun{
#' # build and normalize adjacency
#' adj <- build_network_adj(edgelist)
#' adj_norm <- normalize_network_adj(adj)
#'
#' # synthetic mutmat
#' mutmat <- matrix(0, nrow = nrow(adj_norm), ncol = 2,
#'                  dimnames = list(rownames(adj_norm), c("P1","P2")))
#' mutmat["TP53","P1"] <- 1
#'
#' # run RWR
#' scores <- rwr_multi(adj_norm, mutmat, alpha = 0.7, verbose = FALSE)
#' }
#' @export
#' @importFrom Matrix rowSums
#' @importFrom progress progress_bar
rwr_multi <- function(network.adj_norm, mutmat,
                      alpha    = 0.7,
                      maxiters = 20000,
                      eps      = 1e-14,
                      verbose  = TRUE) {
  # align mutation matrix to network genes
  genes <- rownames(network.adj_norm)
  mm   <- as.data.frame(mutmat)[match(genes, rownames(mutmat)), , drop = FALSE]
  mm[is.na(mm)] <- 0
  mm <- cbind(GENE = genes, mm)

  patientIDs <- colnames(mm)[-1L]
  nP <- length(patientIDs)

  if (verbose) {
    if (requireNamespace("progress", quietly = TRUE)) {
      pb <- progress::progress_bar$new(
        format = "  RWR [:bar] :percent ETA: :eta",
        total  = nP, clear = FALSE, width = 60
      )
    } else {
      message("Note: install.packages('progress') for a progress bar.")
      verbose <- FALSE
    }
  }

  out <- matrix(0,
                nrow = length(genes),
                ncol = nP,
                dimnames = list(genes, patientIDs))

  for (i in seq_len(nP)) {
    pid <- patientIDs[i]
    p0  <- mm[[pid]]
    p   <- p0

    for (t in seq_len(maxiters)) {
      p_old <- p
      p     <- (1 - alpha) * as.vector(network.adj_norm %*% p_old) +
        alpha * p0
      if (sum(abs(p - p_old)) < eps) break
    }

    out[, i] <- p
    if (verbose) {
      if (t == maxiters) {
        message(
          sprintf("Reached maxiters (%d) for patient %s without convergence.",
                  maxiters, pid)
        )
      }
      pb$tick()
    }
  }

  out
}
