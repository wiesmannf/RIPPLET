# Portions adapted from PMAPscore::newspia (Han et al., 2022),
# licensed GPL (>= 2). Original copyright © the PMAPscore authors.
# Modifications © 2025 Fabius Wiesmann.

#' Patient-Specific Pathway Scores (PPS)
#'
#' Compute one pathway score per pathway and per patient from gene-level scores
#' and a list of RIPPLET-normalized pathway adjacency matrices. Returns a
#' ranked table for each patient. Accepts either the full output of
#' \code{\link{get_ripplet_pathways}()} or a plain named list of normalized
#' pathway matrices. Function ‘newspia‘ This function is based on SPIA algorithm, adapted from \code{PMAPscore::newspia} (Han et al., 2022).
#'
#' @param gene_scores Numeric matrix (genes x samples). Row names must be gene
#'   symbols; column names are patient IDs.
#' @param pathways Either the list returned by \code{get_ripplet_pathways()}
#'   (uses \code{$SPIA_Norm_Adj_Matrices}) or a named list of normalized
#'   pathway adjacency matrices (genes x genes).
#' @param threshold_level_ps Numeric scalar; values in \code{gene_scores} below
#'   this threshold are set to 0 (default \code{0.01}).
#' @param show_progress Logical; if \code{TRUE}, show a progress bar when the
#'   \pkg{progress} package is available (default \code{TRUE}).
#'
#' @return A named list with one \code{data.frame} per sample. Each data frame
#'   has columns \code{Pathway} and \code{sumpfs} (sum of perturbation factors), sorted decreasing by \code{sumpfs}.
#'
#' @examples
#' \dontrun{
#' # Gene scores (e.g., from ripplet_gene)
#' # gis <- ripplet_gene(edgelist, mutmat, verbose = FALSE)
#'
#' # Pathways
#' # pw <- get_ripplet_pathways(data.base = "reactome", show_progress = FALSE)
#'
#' # Sample-specific pathway scores
#' # pps <- get_pps(gene_scores = gis, pathways = pw, threshold_level_ps = 0.01)
#' # head(pps[[1]])
#' }
#'
#' @export
get_pps <- function(gene_scores,
                    pathways,
                    threshold_level_ps = 0.01,
                    show_progress = TRUE) {
  # Normalize inputs -----------------------------------------------------------
  if (is.data.frame(gene_scores)) gene_scores <- as.matrix(gene_scores)
  if (!is.matrix(gene_scores))
    stop("'gene_scores' must be a numeric matrix (genes x samples).")
  if (is.null(rownames(gene_scores)))
    stop("'gene_scores' must have row names (gene symbols).")
  if (is.null(colnames(gene_scores)))
    stop("'gene_scores' must have column names (sampleIDs).")

  # Accept either full get_ripplet_pathways() output or a plain list of matrices
  if (is.list(pathways) && "SPIA_Norm_Adj_Matrices" %in% names(pathways)) {
    adj_list <- pathways$SPIA_Norm_Adj_Matrices
  } else {
    adj_list <- pathways
  }
  if (!length(adj_list))
    stop("'pathways' (or $SPIA_Norm_Adj_Matrices) appears to be empty.")

  # Threshold gene scores ------------------------------------------------------
  gs <- gene_scores
  gs[gs < threshold_level_ps] <- 0

  # Output container -----------------------------------------------------------
  n_pat <- ncol(gs)
  pps_res.list <- vector("list", n_pat)
  names(pps_res.list) <- colnames(gs)

  # Optional progress bar ------------------------------------------------------
  pb <- NULL
  if (show_progress) {
    if (requireNamespace("progress", quietly = TRUE)) {
      pb <- progress::progress_bar$new(
        format = "  Processing [:bar] :current/:total (:percent) eta: :eta",
        total  = n_pat, clear = FALSE, width = 60
      )
    } else {
      message("Note: install.packages('progress') for a progress bar.")
      show_progress <- FALSE
    }
  }

  # Main loop over samples ----------------------------------------------------
  for (i in seq_len(n_pat)) {
    if (show_progress) pb$tick()

    sampleID <- colnames(gs)[i]
    scores_i  <- gs[, i]

    # Vector to hold sumpfs per pathway
    nPath <- length(adj_list)
    sumpfs   <- numeric(nPath)
    names(sumpfs) <- names(adj_list)

    # Compute sumpfs for each pathway
    for (j in seq_along(adj_list)) {
      M <- adj_list[[j]]
      diag(M) <- diag(M) - 1L

      X <- scores_i[rownames(M)]
      X[is.na(X)] <- 0

      # Solve; treat non-invertible systems as NA
      pf <- tryCatch(solve(M, -X), error = function(e) NULL)
      if (is.null(pf)) {
        sumpfs[j] <- NA_real_
      } else {
        sumpfs[j] <- sum(pf, na.rm = TRUE)
      }
    }

    # Assemble and sort results
    res <- data.frame(
      Pathway = names(sumpfs),
      sumpfs     = sumpfs,
      stringsAsFactors = FALSE
    )
    res <- res[order(res$sumpfs, decreasing = TRUE), ]
    rownames(res) <- NULL

    pps_res.list[[sampleID]] <- res
  }

  pps_res.list
}


#' Normalized Pathway Scores (NPPS)
#'
#' Convert sample-specific pathway tables (from \code{get_pps()}) into a
#' pathway-by-sample matrix of normalized scores by dividing each observed
#' pathway score by its corresponding MaxHit value.
#'
#' @param pps A named list of data frames as returned by \code{get_pps()}.
#'   Each data frame must contain columns \code{Pathway} and
#'   score column (\code{"sumpfs"}).
#' @param pathways Either the full list returned by
#'   \code{\link{get_ripplet_pathways}()} (uses \code{$MaxHit}), or a named
#'   numeric vector of MaxHit scores (names are pathway IDs).
#' @param show_progress Logical; if \code{TRUE}, show a progress bar when the
#'   \pkg{progress} package is available (default \code{TRUE}).
#'
#' @return A numeric matrix of normalized pathway scores (rows = pathways,
#'   columns = samples). Entries are \code{NA} if the pathway is missing in
#'   a sample or if MaxHit is missing/non-positive.
#'
#' @examples
#' \dontrun{
#' # Gene scores:
#' # gis <- ripplet_gene(edgelist, mutmat, verbose = FALSE)
#' # Pathways:
#' # pw  <- get_ripplet_pathways(data.base = "reactome", show_progress = FALSE)
#' # Per-patient pathway tables:
#' # pps <- get_pps(gene_scores = gis, pathways = pw, threshold_level_ps = 0.01)
#' # Normalized matrix:
#' # npps <- get_npps(pps = pps, pathways = pw)
#' # dim(npps); npps[1:5, 1:3]
#' }
#'
#' @export
get_npps <- function(pps,
                     pathways,
                     show_progress = TRUE) {
  # Validate pps ---------------------------------------------------------------
  if (!is.list(pps) || !length(pps))
    stop("'pps' must be a non-empty named list of data frames (from get_pps()).")
  if (is.null(names(pps)) || any(!nzchar(names(pps))))
    stop("'pps' must be a named list; names should be sampleIDs.")
  # Check that required columns exist in first element (cheap check)
  first_df <- pps[[1L]]
  if (!is.data.frame(first_df) ||
      !"Pathway" %in% names(first_df)) {
    stop("Each element of 'pps' must be a data frame with columns ",
         "'Pathway'")
  }

  # Extract MaxHit vector ------------------------------------------------------
  if (is.list(pathways) && "MaxHit" %in% names(pathways)) {
    maxhit <- pathways$MaxHit
  } else if (is.numeric(pathways) && !is.null(names(pathways))) {
    maxhit <- pathways
  } else {
    stop("'pathways' must be either the list returned by get_ripplet_pathways() ",
         "or a named numeric vector of MaxHit scores.")
  }
  if (!length(maxhit))
    stop("'MaxHit' appears to be empty.")

  # Determine pathway universe (intersection is safest) -----------------------
  pps_paths <- unique(unlist(lapply(pps, function(df) as.character(df$Pathway))))
  all_paths <- intersect(pps_paths, names(maxhit))
  if (!length(all_paths))
    stop("No overlapping pathways between 'pps' and 'MaxHit' names.")
  # Inform if anything is dropped
  dropped_pps <- setdiff(pps_paths, all_paths)
  dropped_mxh <- setdiff(names(maxhit), all_paths)
  if (length(dropped_pps))
    message("Note: dropping ", length(dropped_pps),
            " pps pathways not present in MaxHit.")
  if (length(dropped_mxh))
    message("Note: ", length(dropped_mxh),
            " MaxHit pathways not present in pps; they will not appear in output.")

  # Initialize output matrix ---------------------------------------------------
  samples <- names(pps)
  npps_mat <- matrix(NA_real_,
                     nrow = length(all_paths),
                     ncol = length(samples),
                     dimnames = list(all_paths, samples))

  # Optional progress bar ------------------------------------------------------
  pb <- NULL
  if (show_progress) {
    if (requireNamespace("progress", quietly = TRUE)) {
      pb <- progress::progress_bar$new(
        format = "  Normalizing [:bar] :current/:total (:percent) eta: :eta",
        total  = length(samples), clear = FALSE, width = 60
      )
    } else {
      message("Note: install.packages('progress') for a progress bar.")
      show_progress <- FALSE
    }
  }

  # Fill matrix ----------------------------------------------------------------
  for (pid in samples) {
    if (show_progress) pb$tick()
    df <- pps[[pid]]
    # Keep only pathways we will output
    df <- df[df$Pathway %in% all_paths, , drop = FALSE]

    if (!nrow(df)) next

    obs <- df[["sumpfs"]]
    pth <- as.character(df$Pathway)

    # Align denominators
    denom <- maxhit[pth]
    # Valid where denominator is finite and > 0
    valid <- is.finite(obs) & is.finite(denom) & denom > 0

    norm_vals <- rep(NA_real_, length(pth))
    norm_vals[valid] <- obs[valid] / denom[valid]

    # Place into matrix rows matching pathways
    npps_mat[pth, pid] <- norm_vals
  }

  npps_mat
}
