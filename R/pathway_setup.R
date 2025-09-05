# SPDX-License-Identifier: AGPL-3.0-or-later
# This file is part of RIPPLET (AGPL-3.0-or-later).
# Portions adapted from PMAPscore (GPL-3.0-or-later). Copyright © PMAPscore authors.
# We selected the GPL-3 option for the adapted parts. See inst/COPYRIGHTS for details.


#' Retrieve Pathway Graphs for RIPPLET
#'
#' Download pathway graphs from a selected database via \pkg{graphite}
#' for human ("hsapiens") and convert node identifiers to HGNC symbols.
#'
#' @param data.base Character scalar, pathway database to use.
#'   Options include "kegg", "panther", "pathbank", "pharmgkb", "reactome", "smpdb",
#'   "wikipathways", etc. See \code{graphite::pathwayDatabases()}.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{data.base} - the database name.
#'   \item \code{pways} - a list of \pkg{graphite} pathway objects with
#'         node IDs converted to gene symbols.
#' }
#'
#' @examples
#' \dontrun{
#' pw <- get_graphs("reactome")
#' length(pw$pways)
#' }
#'
#' @export
#' @importFrom graphite pathways convertIdentifiers
get_graphs <- function(
    data.base = "reactome"
) {
  message("1/2 - Loading & converting identifiers for database: ", data.base)
  pways <- graphite::pathways("hsapiens", data.base)
  pways <- graphite::convertIdentifiers(pways, "SYMBOL")
  list(data.base = data.base, pways = pways)
}


#' Normalize Pathway Graphs for RIPPLET
#'
#' Convert a list of pathway graphs (from \code{\link{get_graphs}})
#' into SPIA-style adjacency matrices and apply normalization so that
#' each matrix can be used for RIPPLET pathway propagation.
#'
#' @param pathway_input Output of \code{\link{get_graphs}}.
#' @param dampening_alpha Numeric in (0,1]; dampening factor applied when
#'   a pathway matrix is close to singular (default 0.99).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{INFO} - summary table of processed pathways.
#'   \item \code{Adj_Matrices} - raw signed adjacency matrices.
#'   \item \code{SPIA_Norm_Adj_Matrices} - normalized adjacency matrices.
#' }
#'
#' @examples
#' \dontrun{
#' pw <- get_graphs("reactome")
#' norm <- norm_graphs(pw)
#' head(norm$INFO)
#' }
#'
#' @export
#' @importFrom utils getFromNamespace
norm_graphs <- function(
    pathway_input,
    dampening_alpha = 0.99
) {
  data.base <- pathway_input$data.base
  pways     <- pathway_input$pways

  # fetch graphite internals
  translateEdges <- utils::getFromNamespace("translateEdges", "graphite")
  edgeMatrix     <- utils::getFromNamespace("edgeMatrix",     "graphite")
  spiaAttributes <- utils::getFromNamespace("spiaAttributes", "graphite")

  # Relations & weights
  rel <- c(
    "activation","compound","binding/association","expression","inhibition",
    "activation_phosphorylation","phosphorylation","inhibition_phosphorylation",
    "inhibition_dephosphorylation","dissociation","dephosphorylation",
    "activation_dephosphorylation","state change","activation_indirect effect",
    "inhibition_ubiquination","ubiquination","expression_indirect effect",
    "inhibition_indirect effect","repression","dissociation_phosphorylation",
    "indirect effect_phosphorylation","activation_binding/association",
    "indirect effect","activation_compound","activation_ubiquination"
  )
  beta <- c(1,0,0,1,1,1,0,1,1,0,0,1,0,1,1,0,1,1,1,0,0,1,0,1,1)
  names(beta) <- rel

  message("2/2 - Preparing SPIA & normalizing for database: ", data.base)

  path.info <- lapply(pways, function(pw) {
    tl <- translateEdges(pw)
    if (is.null(tl) || nrow(tl$edges) < 5) return(NULL)
    f <- edgeMatrix(tl$nodes, tl$edges)
    mats <- lapply(spiaAttributes, f)
    names(mats) <- spiaAttributes
    mats$title             <- pw@title
    mats$nodes             <- tl$nodes
    mats$NumberOfReactions <- 0L
    mats
  })
  path.info <- Filter(Negate(is.null), path.info)
  names(path.info) <- vapply(path.info, `[[`, "", "title")
  datpT <- path.info

  nP           <- length(datpT)
  adj          <- vector("list", nP)
  datp         <- vector("list", nP)
  path.names   <- character(nP)
  hasR         <- logical(nP)
  sizem        <- integer(nP)
  damped_flag  <- logical(nP)
  names(adj)   <- names(datp) <- names(datpT)

  for (jj in seq_along(datpT)) {
    s   <- 0; con <- 0
    for (r in rel) {
      con <- con + datpT[[jj]][[r]] * abs(sign(beta[r]))
      s   <- s   + datpT[[jj]][[r]] * beta[r]
    }
    sizem[jj]     <- nrow(s)
    adj[[jj]]     <- s

    z <- matrix(rep(colSums(con), times = nrow(con)),
                nrow = nrow(con), ncol = ncol(con), byrow = TRUE)
    z[z == 0] <- 1

    M <- s / z

    Mtmp <- M; diag(Mtmp) <- diag(Mtmp) - 1
    if (abs(det(Mtmp)) <= 1e-7) {
      M <- dampening_alpha * M
      damped_flag[jj] <- TRUE
    }
    datp[[jj]]     <- M
    path.names[jj] <- datpT[[jj]]$title
    hasR[jj]       <- datpT[[jj]]$NumberOfReactions >= 1
  }

  clean_names <- function(mat) {
    rownames(mat) <- sub("^SYMBOL:", "", rownames(mat))
    colnames(mat) <- sub("^SYMBOL:", "", colnames(mat))
    mat
  }
  adj  <- lapply(adj,  clean_names)
  datp <- lapply(datp, clean_names)

  zero_sum_or_noRx <-
    (vapply(datp, function(m) sum(abs(m)), numeric(1)) == 0) |
    hasR |
    is.na(path.names)
  keep1 <- !zero_sum_or_noRx
  adj_pre         <- adj[keep1]
  datp_pre        <- datp[keep1]
  names_pre       <- path.names[keep1]
  damped_pre      <- damped_flag[keep1]
  sizem_pre       <- sizem[keep1]

  path.invert <- vapply(datp_pre, function(M) {
    M2 <- M; diag(M2) <- diag(M2) - 1
    abs(det(M2)) > 1e-7
  }, logical(1))
  names(path.invert) <- names(datp_pre)
  adj_fin  <- adj_pre[path.invert]
  datp_fin <- datp_pre[path.invert]

  INFO <- data.frame(
    database   = data.base,
    pathway    = names_pre,
    psize      = sizem_pre,
    damped     = damped_pre,
    damp_fact  = ifelse(damped_pre, dampening_alpha, NA),
    invertible = path.invert,
    stringsAsFactors = FALSE
  )
  INFO$final_pass <- INFO$invertible

  list(
    INFO                   = INFO,
    Adj_Matrices           = adj_fin,
    SPIA_Norm_Adj_Matrices = datp_fin
  )
}

#' Maximum Pathway Propagation Score
#'
#' For each normalized pathway adjacency matrix, compute a maximum
#' propagation score across all nodes. Useful as a simple summary of
#' pathway activity. Matrices that are not invertible return \code{NA}.
#'
#' @param in_adj_mat_norm A named list of normalized pathway adjacency matrices.
#' @param show_progress Logical; if \code{TRUE}, show a progress bar
#'   (requires the \pkg{progress} package). Default is \code{TRUE}.
#'
#' @return A named numeric vector of maximum scores per pathway,
#'   sorted in decreasing order.
#'
#' @examples
#' \dontrun{
#' scores <- get_pps_max_hit(norm_list, show_progress = FALSE)
#' head(scores)
#' }
#'
#' @export
get_pps_max_hit <- function(in_adj_mat_norm, show_progress = TRUE) {
  if (show_progress) {
    if (requireNamespace("progress", quietly = TRUE)) {
      total_iter <- sum(vapply(in_adj_mat_norm, nrow, integer(1)))
      pb <- progress::progress_bar$new(
        total  = total_iter,
        format = "Processing [:bar] :percent ETA: :eta"
      )
    } else {
      message("Note: install.packages('progress') for a progress bar.")
      show_progress <- FALSE
    }
  }

  results <- vector("list", length(in_adj_mat_norm))
  names(results) <- names(in_adj_mat_norm)

  for (p in seq_along(in_adj_mat_norm)) {
    M <- in_adj_mat_norm[[p]]
    diag(M) <- diag(M) - 1
    n <- nrow(M)

    if (n == 0L) {
      sumpfs <- numeric(0)
      if (show_progress) pb$tick(0)
    } else {
      # Single factorization: solve once with matrix RHS (columns = unit impulses)
      B  <- -diag(n)
      PF <- tryCatch(solve(M, B), error = function(e) NULL)

      if (is.null(PF)) {
        sumpfs <- rep(NA_real_, n)
        if (show_progress) pb$tick(n)
      } else {
        sumpfs <- colSums(PF, na.rm = TRUE)
        if (show_progress) pb$tick(n)
      }
    }

    maxSumpfs <- if (length(sumpfs) == 0L || all(is.na(sumpfs))) NA_real_ else max(sumpfs, na.rm = TRUE)
    results[[p]] <- list(sumPfs = sumpfs, maxSumpfs = maxSumpfs)
  }

  maxhit_scores <- vapply(results, `[[`, numeric(1), "maxSumpfs")
  sort(maxhit_scores, decreasing = TRUE)
}
