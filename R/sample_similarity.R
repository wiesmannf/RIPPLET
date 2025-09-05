#' Compute Sample Similarity from RWR Profiles
#'
#' Given a matrix or data.frame of Random‑Walk‑with‑Restart steady‑state
#' scores (genes × samples), this function
#' 1. normalizes each column to sum to 1,
#' 2. log‑transforms with a small constant to avoid division by zero,
#' 3. performs a truncated SVD (via \code{irlba}) of rank \code{k},
#'    with the RNG seeded locally for reproducibility,
#' 4. maps the left singular vectors into a low‑dimensional space,
#' 5. computes angular similarity between sample vectors.
#'
#' @param RWR_patspec A numeric matrix or data.frame, rows = genes,
#'   columns = samples (e.g. output from \code{\link{rwr_multi}()}).
#' @param k Integer; number of singular vectors to retain (default 20).
#' @param seed Integer; seed for the truncated SVD RNG (default 1).
#' @return A symmetric numeric matrix of angular similarities
#'   (samples × samples), with row‐ and column‐names from \code{colnames(RWR_patspec)}.
#' @examples
#' \dontrun{
#' # Suppose `scores` is your RWR_patspec matrix:
#' sim <- compute_sample_similarity(scores, k = 10, seed = 42)
#' }
#' @export
#' @importFrom withr with_seed
#' @importFrom irlba irlba
#' @importFrom lsa cosine
compute_sample_similarity <- function(RWR_patspec, k = 20, seed = 1) {
  # Convert to matrix and normalize columns
  mat <- as.matrix(RWR_patspec)
  col_sums <- colSums(mat, na.rm = TRUE)
  col_sums[col_sums == 0] <- 1
  mat_norm <- sweep(mat, 2, col_sums, FUN = "/")

  # Log transform with small constant
  F_mat <- t(mat_norm)
  c_const <- 1 / ncol(F_mat)
  F_log <- log(F_mat + c_const)

  # Safe k
  k_max <- max(1L, min(nrow(F_log), ncol(F_log)) - 1L)
  if (k > k_max) {
    message("Requested k = ", k,
            " is too large; using maximum allowable k = ", k_max, " instead.")
    k <- k_max
  }

  # Truncated SVD with reproducible seed
  svd_res <- withr::with_seed(seed,
                              irlba::irlba(F_log, nv = k)
  )
  U    <- svd_res$u
  sing <- svd_res$d

  # Build low‑dim representation
  S_half <- diag(sqrt(sing), nrow = k, ncol = k)
  M      <- S_half %*% t(U)

  # Cosine -> angular similarity
  cos_sim      <- lsa::cosine(M)
  angular_dist <- acos(cos_sim)
  angular_sim  <- 1 - (angular_dist / pi)

  # Preserve sample names
  sample_names      <- colnames(mat_norm)
  rownames(angular_sim) <- sample_names
  colnames(angular_sim) <- sample_names

  angular_sim
}
