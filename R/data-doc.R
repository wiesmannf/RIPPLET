#' Example edge list of a small PPI network
#'
#' Undirected edges between gene symbols used in examples/tests.
#'
#' @format A data frame with 3 columns:
#' \describe{
#'   \item{From}{Source gene symbol (character).}
#'   \item{To}{Target gene symbol (character).}
#'   \item{weight}{Edge weight (numeric, usually 1).}
#' }
#' @usage data(edgelist)
#' @keywords datasets
#' @examples
#' data(edgelist)
#' head(edgelist)
"edgelist"


#' Example binary mutation matrix
#'
#' Gene (rows) × patient (columns) 0/1 mutation matrix for examples/tests.
#'
#' @format A matrix with genes as rownames and patient IDs as colnames. Entries are 0/1.
#' @usage data(mutmat)
#' @keywords datasets
#' @examples
#' data(mutmat)
#' dim(mutmat)
"mutmat"

#' Reference gene set used in tests
#'
#' Character vector of gene symbols used as a reference set.
#'
#' @format A character vector of gene symbols.
#' @usage data(refgenes)
#' @keywords datasets
#' @examples
#' data(refgenes)
#' length(refgenes)
"refgenes"
