#' Compute Feature×Patient Scores Weighted by Patient Similarity
#'
#' Given a matrix of per‐feature RWR scores (\code{genes × patients}) and
#' a square patient–patient similarity matrix (default: angular similarity), this function for each patient
#' computes a weighted aggregate of other patients’ profiles (using a Gaussian‐
#' like decay on similarity), multiplies it elementwise by the index patient’s
#' own profile, and rescales each result to \[0,1\]. This results in the gene impact matrix.
#'
#' @param RWR_patspec Numeric matrix (features × patients) of RWR scores.
#' @param sim_mat Numeric square matrix (patients × patients) of similarities
#'   in \[0,1\], with identical row and column names.
#' @param delta Numeric scalar > 0; controls the width of the Gaussian‐like
#'   weight decay (default 0.2).
#' @return A numeric matrix of the same dimensions and dimnames as
#'   \code{RWR_patspec}, containing gene impact scores.
#' @examples
#' \dontrun{
#' # simple RWR profiles for 2 patients over 3 features
#' patspec <- matrix(c(1,0,0, 0,1,0), nrow=3,
#'                   dimnames=list(paste0("F",1:3), c("P1","P2")))
#' # identity similarity
#' sim <- diag(1,2); rownames(sim) <- colnames(sim) <- c("P1","P2")
#' impact_scoring(patspec, sim)
#' }
#' @export
impact_scoring <- function(RWR_patspec, sim_mat, delta = 0.2) {
  # 0. Sanity checks on sim_mat
  if (!is.matrix(sim_mat) || nrow(sim_mat) != ncol(sim_mat)) {
    stop("sim_mat must be a square matrix")
  }
  if (is.null(rownames(sim_mat)) || is.null(colnames(sim_mat))) {
    stop("sim_mat must have both row names and column names as sample names")
  }
  if (!all(rownames(sim_mat) == colnames(sim_mat))) {
    stop("Row names and column names of sim_mat must be identical and in the same order")
  }

  # Align sim_mat to patients in RWR_patspec
  pats   <- colnames(RWR_patspec)
  sim_ids <- rownames(sim_mat)
  missing <- setdiff(pats, sim_ids)
  if (length(missing) > 0) {
    stop(sprintf(
      "sim_mat is missing similarity data for patients: %s",
      paste(missing, collapse = ", ")
    ))
  }
  extras <- setdiff(sim_ids, pats)
  if (length(extras) > 0) {
    message(sprintf(
      "Extra patients detected in sim_mat; excluded: %s",
      paste(extras, collapse = ", ")
    ))
  }
  sim_mat <- sim_mat[pats, pats]

  # 1. Prepare output container
  n_pat  <- length(pats)
  n_feat <- nrow(RWR_patspec)
  out    <- matrix(
    0,
    nrow    = n_feat,
    ncol    = n_pat,
    dimnames= list(rownames(RWR_patspec), pats)
  )

  # Weighting function
  weight_fn <- function(sim) {
    exp(- (1 - sim)^2 / (2 * delta^2))
  }

  # 2. Loop over patients
  for (i in seq_along(pats)) {
    sim_vec <- sim_mat[, i]
    sim_vec[i] <- 0                   # self-similarity to 0
    w <- weight_fn(sim_vec)
    w <- w / sum(w)                   # normalize

    agg_pattern <- as.vector(RWR_patspec %*% w)
    raw_score   <- RWR_patspec[, i] * agg_pattern

    rng <- range(raw_score, na.rm = TRUE)
    if (diff(rng) == 0) {
      scaled <- rep(0, length(raw_score))
    } else {
      scaled <- (raw_score - rng[1]) / diff(rng)
    }

    out[, i] <- scaled
  }

  out
}
