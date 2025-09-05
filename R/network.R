#' Build a Sparse (Undirected) Adjacency Matrix from an Edge List
#'
#' Creates an undirected igraph graph from a two-column edge list,
#' removes loops and multiple edges, and returns its \code{Matrix::dgCMatrix}
#' adjacency matrix (symmetrical, unweighted).
#'
#' @param edgelist A \code{data.frame} or \code{tibble} with two columns
#'   (source, target), each identifying a node.
#' @return A sparse adjacency matrix of class \linkS4class{dgCMatrix}.
#' @examples
#' \dontrun{
#' df <- data.frame(from = c("A", "A", "B"), to = c("B", "B", "C"))
#' adj <- build_network_adj(df)
#' }
#' @export
#' @importFrom igraph graph_from_data_frame simplify as_adj
build_network_adj <- function(edgelist) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Package 'igraph' is required but not installed.", call. = FALSE)

  g <- igraph::graph_from_data_frame(edgelist, directed = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

  # Return a sparse dgCMatrix with vertex names as dimnames (same as before)
  igraph::as_adj(g, names = TRUE, sparse = TRUE, type = "both")
}



#' Row-Normalize a Sparse Adjacency Matrix
#'
#' Scales each row of a sparse adjacency matrix so that non-empty rows sum to 1,
#' producing a row-stochastic matrix while leaving empty rows as all-zeros. This
#' matrix is to be used for the Random Walk with Restart (RWR).
#'
#' @param network.adj A sparse adjacency matrix (e.g. a
#'   \linkS4class{dgCMatrix}) with nonnegative entries.
#' @return A sparse matrix of the same size and dimnames, row-stochastic
#'   where possible and zero where input rows were empty.
#' @examples
#' \dontrun{
#' adj <- Matrix::sparseMatrix(
#'   i = c(1,1,2,2,3), j = c(2,3,1,3,1), x = 1,
#'   dimnames = list(c("A","B","C"), c("A","B","C"))
#' )
#' adj_norm <- normalize_network_adj(adj)
#'
#' # Each non-empty row now sums to 1
#' rowSums(as.matrix(adj_norm))
#' }
#' @export
#' @importFrom Matrix rowSums Diagonal
normalize_network_adj <- function(network.adj) {
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' is required but not installed.", call. = FALSE)

  # Compute inverse row sums (with safe handling of zero‐sum rows)
  row_sums    <- Matrix::rowSums(network.adj)
  inv_row_sums <- ifelse(row_sums > 0, 1 / row_sums, 0)

  # Build diagonal matrix of inverse row sums and multiply
  D         <- Matrix::Diagonal(x = inv_row_sums)
  adj_norm  <- D %*% network.adj

  # Preserve dimnames
  dimnames(adj_norm) <- dimnames(network.adj)
  adj_norm
}
