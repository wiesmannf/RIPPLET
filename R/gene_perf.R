#' Personalized Ranking Performance from Impact Scores
#'
#' Compute per-sample precision/recall/F1 curves and normalized AUCs (npAUC)
#' using a matrix of gene impact scores (e.g., output of \code{\link{impact_scoring}}),
#' a binary mutation matrix, and a global reference gene set.
#'
#' @param global_ref Character vector of reference genes (e.g., disease/driver set).
#'   See dataset \code{refgenes} for an example.
#' @param mut_mat Binary matrix (genes × patients) of mutations (0/1).
#' @param impact_scores Numeric matrix (genes × patients), typically from
#'   \code{\link{impact_scoring}}; higher scores rank earlier.
#' @param max_rank Optional integer. If \code{NULL}, chosen adaptively from the data.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{avg_df}: data.frame with columns \code{Rank}, \code{Precision}, \code{Recall}, \code{F1}
#'         averaged across patients up to \code{max_rank}.
#'   \item \code{auc_stats}: data.frame with per-metric mean and SD of patient-level npAUCs.
#'   \item \code{max_rank}: the rank cutoff used.
#' }
#'
#' @examples
#' \dontrun{
#' data(refgenes)
#' data(mutmat)
#' # Suppose 'gis' is the impact score matrix from impact_scoring(...)
#' # perf <- personalized_perf(global_ref = refgenes, mut_mat = mutmat, impact_scores = gis)
#' # head(perf$avg_df); perf$auc_stats
#' }
#' @export
personalized_perf <- function(global_ref, mut_mat, impact_scores, max_rank = NULL) {
  # tiny local trapezoid integrator (no external dependency)
  trapz_local <- function(x, y) {
    if (length(x) < 2L) return(0)
    dx <- diff(x)
    # same length as dx
    ym <- (y[-length(y)] + y[-1L]) / 2
    sum(dx * ym, na.rm = TRUE)
  }

  # --- Basic validation -------------------------------------------------------
  if (is.data.frame(mut_mat))        mut_mat        <- as.matrix(mut_mat)
  if (is.data.frame(impact_scores))  impact_scores  <- as.matrix(impact_scores)
  if (!is.matrix(mut_mat) || is.null(rownames(mut_mat)) || is.null(colnames(mut_mat))) {
    stop("'mut_mat' must be a matrix with gene rownames and patient colnames.")
  }
  if (!is.matrix(impact_scores) || is.null(rownames(impact_scores)) || is.null(colnames(impact_scores))) {
    stop("'impact_scores' must be a matrix with gene rownames and patient colnames.")
  }
  if (!is.character(global_ref) || length(global_ref) == 0) {
    stop("'global_ref' must be a non-empty character vector of gene symbols.")
  }

  # --- Align genes & patients -------------------------------------------------
  genes <- intersect(rownames(impact_scores), rownames(mut_mat))
  if (length(genes) == 0) {
    mr <- if (!is.null(max_rank) && is.numeric(max_rank) && length(max_rank) == 1 && max_rank >= 1) {
      as.integer(max_rank)
    } else 20L
    avg_df <- data.frame(
      Rank      = seq_len(mr),
      Precision = rep(NA_real_, mr),
      Recall    = rep(NA_real_, mr),
      F1        = rep(NA_real_, mr),
      stringsAsFactors = FALSE
    )
    auc_stats <- data.frame(
      Metric  = c("Precision", "Recall", "F1"),
      MeanAUC = NA_real_,
      SDAUC   = NA_real_,
      stringsAsFactors = FALSE
    )
    return(list(avg_df = avg_df, auc_stats = auc_stats, max_rank = mr))
  }
  mut_mat       <- mut_mat[genes, , drop = FALSE]
  impact_scores <- impact_scores[genes, , drop = FALSE]

  pats <- intersect(colnames(impact_scores), colnames(mut_mat))
  if (length(pats) == 0) {
    mr <- if (!is.null(max_rank) && is.numeric(max_rank) && length(max_rank) == 1 && max_rank >= 1) {
      as.integer(max_rank)
    } else 20L
    avg_df <- data.frame(
      Rank      = seq_len(mr),
      Precision = rep(NA_real_, mr),
      Recall    = rep(NA_real_, mr),
      F1        = rep(NA_real_, mr),
      stringsAsFactors = FALSE
    )
    auc_stats <- data.frame(
      Metric  = c("Precision", "Recall", "F1"),
      MeanAUC = NA_real_,
      SDAUC   = NA_real_,
      stringsAsFactors = FALSE
    )
    return(list(avg_df = avg_df, auc_stats = auc_stats, max_rank = mr))
  }

  impact_scores <- impact_scores[, pats, drop = FALSE]
  mut_mat       <- mut_mat[,       pats, drop = FALSE]
  num_pat <- length(pats)

  # --- Choose max_rank if not provided ---------------------------------------
  if (is.null(max_rank)) {
    pers_sizes <- vapply(pats, function(p) {
      sum(rownames(mut_mat)[mut_mat[, p] > 0] %in% global_ref)
    }, integer(1))
    valid <- pers_sizes >= 3
    if (!any(valid)) {
      max_rank <- 20L
    } else {
      max_rank <- max(1L, as.integer(2 * stats::median(pers_sizes[valid])))
    }
  } else {
    if (!is.numeric(max_rank) || length(max_rank) != 1 || max_rank < 1) {
      stop("'max_rank' must be a single positive integer or NULL.")
    }
    max_rank <- as.integer(max_rank)
  }

  # --- Storage ---------------------------------------------------------------
  precision_mat <- matrix(NA_real_, nrow = max_rank, ncol = num_pat)
  recall_mat    <- matrix(NA_real_, nrow = max_rank, ncol = num_pat)
  f1_mat        <- matrix(NA_real_, nrow = max_rank, ncol = num_pat)
  npAUC_prec    <- numeric(num_pat)
  npAUC_rec     <- numeric(num_pat)
  npAUC_f1      <- numeric(num_pat)
  names(npAUC_prec) <- pats

  # --- Per-patient loop -------------------------------------------------------
  for (i in seq_along(pats)) {
    p      <- pats[i]
    scores <- impact_scores[, p]

    # Ranking (high score first); break ties by gene name for determinism
    ord     <- order(-scores, rownames(impact_scores))
    ranked  <- rownames(impact_scores)[ord]
    n_ranked <- length(ranked)

    if (n_ranked == 0) {
      npAUC_prec[i] <- NA; npAUC_rec[i] <- NA; npAUC_f1[i] <- NA
      next
    }
    mutated  <- rownames(mut_mat)[mut_mat[, p] > 0]
    pers_ref <- intersect(mutated, global_ref)
    if (length(pers_ref) == 0) {
      npAUC_prec[i] <- NA; npAUC_rec[i] <- NA; npAUC_f1[i] <- NA
      next
    }

    K <- min(n_ranked, max_rank)
    for (k in seq_len(K)) {
      topk <- ranked[1:k]
      TP   <- sum(topk %in% pers_ref)
      FP   <- k - TP
      FN   <- length(pers_ref) - TP
      prec <- if ((TP + FP) > 0) TP / (TP + FP) else 0
      rec  <- if ((TP + FN) > 0) TP / (TP + FN) else 0
      f1   <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0

      precision_mat[k, i] <- prec
      recall_mat[k,    i] <- rec
      f1_mat[k,        i] <- f1
    }

    # Normalized AUC by global max_rank
    ranks <- seq_len(K)
    npAUC_prec[i] <- trapz_local(ranks, precision_mat[1:K, i]) / max_rank
    npAUC_rec[i]  <- trapz_local(ranks, recall_mat[1:K,    i]) / max_rank
    npAUC_f1[i]   <- trapz_local(ranks, f1_mat[1:K,        i]) / max_rank
  }

  # --- Aggregate across patients ---------------------------------------------
  avg_df <- data.frame(
    Rank      = seq_len(max_rank),
    Precision = rowMeans(precision_mat, na.rm = TRUE),
    Recall    = rowMeans(recall_mat,    na.rm = TRUE),
    F1        = rowMeans(f1_mat,        na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  auc_stats <- data.frame(
    Metric  = c("Precision", "Recall", "F1"),
    Mean_npAUC = c(mean(npAUC_prec, na.rm = TRUE),
                mean(npAUC_rec,  na.rm = TRUE),
                mean(npAUC_f1,   na.rm = TRUE)),
    SD_npAUC   = c(stats::sd(npAUC_prec, na.rm = TRUE),
                stats::sd(npAUC_rec,  na.rm = TRUE),
                stats::sd(npAUC_f1,   na.rm = TRUE)),
    stringsAsFactors = FALSE
  )

  list(
    avg_df    = avg_df,
    auc_stats = auc_stats,
    max_rank  = max_rank
  )
}
