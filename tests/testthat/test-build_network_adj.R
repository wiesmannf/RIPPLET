# tests/testthat/test-build_network_adj.R

testthat::test_that("returns a dgCMatrix for a small undirected graph", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")

  df <- data.frame(from = c("A","A","B"), to = c("B","B","C"), stringsAsFactors = FALSE)
  adj <- build_network_adj(df)

  testthat::expect_s4_class(adj, "dgCMatrix")
})

testthat::test_that("matrix is symmetric (undirected)", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")

  df <- data.frame(from = c("A","B"), to = c("B","C"), stringsAsFactors = FALSE)
  adj <- build_network_adj(df)

  testthat::expect_true(Matrix::isSymmetric(adj))
})

testthat::test_that("loops are removed (zero diagonal)", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")

  df <- data.frame(from = c("A","B","C"), to = c("B","C","C"), stringsAsFactors = FALSE) # C-C loop
  adj <- build_network_adj(df)

  # Coerce to dense for a tiny matrix to safely inspect the diagonal.
  m <- as.matrix(adj)
  testthat::expect_true(all(diag(m) == 0))
})

testthat::test_that("duplicate edges are collapsed (A-B once)", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")

  df <- data.frame(from = c("A","A"), to = c("B","B"), stringsAsFactors = FALSE) # duplicate A-B
  adj <- build_network_adj(df)

  m <- as.matrix(adj)
  testthat::expect_equal(m["A","B"], 1)
})

testthat::test_that("dimnames are node names and match in rows/cols", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")

  df <- data.frame(from = c("A","B"), to = c("B","C"), stringsAsFactors = FALSE)
  adj <- build_network_adj(df)

  testthat::expect_equal(rownames(adj), colnames(adj))
  testthat::expect_equal(sort(rownames(adj)), c("A","B","C"))
})

testthat::test_that("empty edgelist returns a 0x0 sparse matrix", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")

  empty_df <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)
  adj <- build_network_adj(empty_df)

  testthat::expect_s4_class(adj, "dgCMatrix")
  testthat::expect_equal(dim(adj), c(0L, 0L))
})

testthat::test_that("tiny graph has expected adjacency pattern", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")

  df <- data.frame(from = c("A","B","C"), to = c("B","C","D"), stringsAsFactors = FALSE)
  adj <- build_network_adj(df)

  expected_dense <- matrix(
    c(0,1,0,0,
      1,0,1,0,
      0,1,0,1,
      0,0,1,0),
    nrow = 4, byrow = TRUE,
    dimnames = list(c("A","B","C","D"), c("A","B","C","D"))
  )

  testthat::expect_equal(as.matrix(adj)[rownames(expected_dense), colnames(expected_dense)],
                         expected_dense)
})

testthat::test_that("invalid input (not two columns) errors", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")

  bad_df <- data.frame(a = 1:3)
  testthat::expect_error(build_network_adj(bad_df))
})


#-------------------------------------------------------------------------------
### Tests using edgelist data from data slot

testthat::test_that("build_network_adj returns a sparse dgCMatrix from packaged edgelist", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")

  data("edgelist")  # dataset must be in your package's data/ directory
  adj <- build_network_adj(edgelist)

  testthat::expect_s4_class(adj, "dgCMatrix")
})

testthat::test_that("adjacency is symmetric (undirected graph)", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")

  data("edgelist")
  adj <- build_network_adj(edgelist)

  testthat::expect_true(Matrix::isSymmetric(adj))
})


testthat::test_that("diagonal is zero (loops removed)", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")

  data("edgelist")
  adj <- build_network_adj(edgelist)

  n <- nrow(adj)
  # test diagonal without coercing to dense:
  testthat::expect_equal(sum(adj * Matrix::Diagonal(n)), 0)
})

testthat::test_that("row/col names equal the set of nodes in edgelist", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")

  data("edgelist")
  adj <- build_network_adj(edgelist)

  cols <- names(edgelist)[1:2]
  nodes <- sort(unique(unlist(edgelist[, cols, drop = FALSE])))

  # same names in rows & cols; order doesn't have to be sorted
  testthat::expect_equal(rownames(adj), colnames(adj))
  testthat::expect_setequal(rownames(adj), nodes)
})

testthat::test_that("duplicate edges are collapsed and weights are binary", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")

  data("edgelist")
  adj <- build_network_adj(edgelist)

  # 1) all nonzeros should be 1 after igraph::simplify(remove.multiple=TRUE)
  testthat::expect_true(length(adj@x) == 0 || all(adj@x == 1))

  # 2) nonzero count corresponds to the number of unique undirected, non-loop edges
  cols <- names(edgelist)[1:2]
  pairs <- edgelist[, cols, drop = FALSE]
  names(pairs) <- c("a", "b")

  # remove loops
  pairs <- pairs[pairs$a != pairs$b, , drop = FALSE]

  # canonicalize unordered pairs (a,b) == (b,a)
  canon <- ifelse(pairs$a <= pairs$b,
                  paste(pairs$a, pairs$b, sep = "||"),
                  paste(pairs$b, pairs$a, sep = "||"))
  expected_undirected_edges <- length(unique(canon))

  # In an undirected adjacency, each edge appears twice (i,j) and (j,i)
  testthat::expect_equal(Matrix::nnzero(adj) / 2, expected_undirected_edges)
})
