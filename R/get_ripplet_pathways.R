#' RIPPLET Pathway Pipeline
#'
#' Convenience wrapper that:
#' (1) fetches pathways via \pkg{graphite},
#' (2) creates (normalized) adjacency matrices for topological pathway scoring, and
#' (3) computes per-pathway Maximum Propagation Scores (MaxHit) for downstream pathway perturbation score normalization (PPS to NPPS)
#'
#' @param data.base Character scalar, pathway database (e.g., "reactome", "kegg"). See \code{graphite::pathwayDatabases()} to see all databases.
#'   Passed to \code{\link{get_graphs}}.
#' @param dampening_alpha Numeric in (0,1]; dampening factor for near-singular
#'   matrices in \code{\link{norm_graphs}} (default 0.99).
#' @param show_progress Logical; if \code{TRUE}, show a progress bar for MaxHit
#'   (requires the \pkg{progress} package). Default \code{TRUE}.
#'
#' @return A list in the same format as \code{\link{norm_graphs}} with
#'   an additional final element:
#'   \itemize{
#'     \item \code{INFO}
#'     \item \code{Adj_Matrices}
#'     \item \code{SPIA_Norm_Adj_Matrices}
#'     \item \code{MaxHit} (named numeric vector of maximum propagation scores)
#'   }
#'
#' @examples
#' \dontrun{
#' out <- get_ripplet_pathways(data.base = "reactome",
#'                           dampening_alpha = 0.99,
#'                           show_progress = FALSE)
#' head(out$INFO)
#' head(out$MaxHit)
#' }
#'
#' @export
get_ripplet_pathways <- function(data.base = "reactome",
                               dampening_alpha = 0.99,
                               show_progress = TRUE) {
  # 1) Get pathways
  pw <- get_graphs(data.base = data.base)

  # 2) Normalize (SPIA-style)
  norm <- norm_graphs(pathway_input = pw,
                              dampening_alpha = dampening_alpha)

  # 3) MaxHit scores on normalized matrices
  maxhit <- get_pps_max_hit(norm$SPIA_Norm_Adj_Matrices,
                         show_progress = show_progress)

  # 4) Return same structure as norm_graphs + MaxHit (as last item)
  norm$MaxHit <- maxhit
  norm
}
