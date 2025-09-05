#after normalization, each non-zero row sums to 1 and zero rows remain 0, dimensions/names are preserved, values are finite and non-negative.
testthat::test_that("normalize_network_adj: row-stochastic on package data adjacency (sampled)", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")
  testthat::skip_on_cran()

  data("edgelist", envir = environment())

  set.seed(11)
  n <- min(5000L, nrow(edgelist))
  samp <- edgelist[sample.int(nrow(edgelist), n), 1:2, drop = FALSE]

  adj <- build_network_adj(samp)
  adj_norm <- normalize_network_adj(adj)

  deg <- Matrix::rowSums(adj)
  rs  <- Matrix::rowSums(adj_norm)

  # Round to 3 decimals per your spec
  rs_rounded <- round(rs, 3)
  testthat::expect_true(all(rs_rounded[deg > 0] == 1))
  testthat::expect_true(all(rs_rounded[deg == 0] == 0))

  # S4 + inheritance-friendly class check
  testthat::expect_s4_class(adj_norm, "sparseMatrix")
  testthat::expect_equal(dim(adj_norm), dim(adj))
  testthat::expect_equal(rownames(adj_norm), rownames(adj))
  testthat::expect_equal(colnames(adj_norm), colnames(adj))

  testthat::expect_true(all(adj_norm@x >= 0))
  testthat::expect_true(all(is.finite(adj_norm@x)))
})

#an artificially added isolated node stays a zero row in the normalized matrix, while all other rows sum to 1
testthat::test_that("normalize_network_adj: zero-degree rows remain all-zero after normalization", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")
  testthat::skip_on_cran()

  data("edgelist", envir = environment())

  # Build adjacency from a sample, then append an isolated node
  set.seed(22)
  n <- min(3000L, nrow(edgelist))
  samp <- edgelist[sample.int(nrow(edgelist), n), 1:2, drop = FALSE]
  adj  <- build_network_adj(samp)

  # Append 1x1 zero block (isolated node) to create a guaranteed zero-degree row
  iso_block <- Matrix::Matrix(0, 1, 1, sparse = TRUE)
  adj_iso <- Matrix::bdiag(adj, iso_block)
  rn <- c(rownames(adj), "ISOLATE_NODE_X")
  dimnames(adj_iso) <- list(rn, rn)

  adj_norm <- normalize_network_adj(adj_iso)

  # Last row corresponds to the isolate: should be all zeros
  last_row_sum <- Matrix::rowSums(adj_norm)[length(rn)]
  testthat::expect_equal(unname(last_row_sum), 0)  # drop the name attribute

  # Non-isolated rows should sum to ~1
  deg <- Matrix::rowSums(adj_iso)
  rs  <- Matrix::rowSums(adj_norm)
  testthat::expect_true(all(round(rs[deg > 0], 3) == 1))
})

#normalization only rescales existing entries and never introduces new nonzeros, so the adjacency’s structure (i, p slots) is identical
testthat::test_that("normalize_network_adj: sparsity pattern (nonzero locations) unchanged", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")
  testthat::skip_on_cran()

  data("edgelist", envir = environment())

  set.seed(33)
  n <- min(4000L, nrow(edgelist))
  samp <- edgelist[sample.int(nrow(edgelist), n), 1:2, drop = FALSE]

  adj <- build_network_adj(samp)
  adj_norm <- normalize_network_adj(adj)

  # Row scaling should not change the nonzero pattern in a sparse matrix
  testthat::expect_identical(adj@i, adj_norm@i)
  testthat::expect_identical(adj@p, adj_norm@p)

  # Values must change for at least one row with degree > 1 (unless trivial)
  if (any(Matrix::rowSums(adj) > 1)) {
    testthat::expect_true(!identical(adj@x, adj_norm@x))
  }
})

# a matrix with no edges normalizes to the same all-zero sparse matrix, with the same dimnames and no nonzeros
testthat::test_that("normalize_network_adj: all-zero matrix stays all-zero", {
  testthat::skip_if_not_installed("Matrix")

  zero_adj <- Matrix::Matrix(0, 3, 3, sparse = TRUE)
  dimnames(zero_adj) <- list(c("A","B","C"), c("A","B","C"))

  z_norm <- normalize_network_adj(zero_adj)

  # Accept any sparseMatrix subclass (dgCMatrix, ddiMatrix, etc.)
  testthat::expect_s4_class(z_norm, "sparseMatrix")

  testthat::expect_equal(dim(z_norm), c(3L, 3L))
  testthat::expect_identical(rownames(z_norm), c("A","B","C"))
  testthat::expect_identical(colnames(z_norm), c("A","B","C"))

  # All entries remain zero
  testthat::expect_equal(Matrix::rowSums(z_norm), rep(0, 3))
  testthat::expect_equal(Matrix::nnzero(z_norm), 0L)
})

# row and column names of the adjacency remain unchanged after normalization
testthat::test_that("normalize_network_adj: dimnames preserved exactly", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("Matrix")
  testthat::skip_on_cran()

  data("edgelist", envir = environment())

  set.seed(44)
  n <- min(2000L, nrow(edgelist))
  samp <- edgelist[sample.int(nrow(edgelist), n), 1:2, drop = FALSE]

  adj <- build_network_adj(samp)
  rn_before <- rownames(adj); cn_before <- colnames(adj)

  adj_norm <- normalize_network_adj(adj)

  testthat::expect_identical(rownames(adj_norm), rn_before)
  testthat::expect_identical(colnames(adj_norm), cn_before)
})
