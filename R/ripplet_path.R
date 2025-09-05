#' RIPPLET Pathway Wrapper
#'
#' Compute sample-specific pathway tables (PPS) and return a normalized
#' pathway-by-sample matrix (NPPS) in one call.
#'
#' @param gene_scores Numeric matrix (genes x samples). Row names are gene
#'   symbols; column names are sampleIDs.
#' @param pathways The full list returned by \code{\link{get_ripplet_pathways}()},
#'   which must include \code{$SPIA_Norm_Adj_Matrices} and \code{$MaxHit}.
#' @param threshold_level_ps Numeric scalar; values in \code{gene_scores} below
#'   this threshold are set to 0 (default \code{0.01}).
#' @param show_progress Logical; if \code{TRUE}, show progress bars when the
#'   \pkg{progress} package is available (default \code{TRUE}).
#'
#' @return A numeric matrix of normalized pathway scores (rows = pathways,
#'   columns = samples).
#'
#' @seealso \code{\link{get_pps}}, \code{\link{get_npps}},
#'   \code{\link{get_ripplet_pathways}}
#'
#' @examples
#' \dontrun{
#' # Gene scores (e.g., from ripplet_gene)
#' # gis <- ripplet_gene(edgelist, mutmat, verbose = FALSE)
#' # Pathways (includes $SPIA_Norm_Adj_Matrices and $MaxHit)
#' # pw  <- get_ripplet_pathways(data.base = "reactome", show_progress = FALSE)
#' # Normalized pathway-by-sample matrix
#' # npps <- ripplet_path(gene_scores = gis, pathways = pw)
#' # dim(npps); npps[1:5, 1:3]
#' }
#'
#' @export
ripplet_path <- function(gene_scores,
                         pathways,
                         threshold_level_ps = 0.01,
                         show_progress = TRUE) {
  pps  <- get_pps(
    gene_scores        = gene_scores,
    pathways           = pathways,
    threshold_level_ps = threshold_level_ps,
    show_progress      = show_progress
  )

  npps <- get_npps(
    pps             = pps,
    pathways        = pathways,
    show_progress   = show_progress
  )

  npps
}
