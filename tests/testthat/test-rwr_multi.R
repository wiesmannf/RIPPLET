# tests/testthat/test-rwr_multi.R

# Verifies the output is a genes × patients matrix with correct dimensions and matching row/column names.
testthat::test_that("rwr_multi returns genes × patients with correct dimnames", {
  testthat::skip_if_not_installed("Matrix")

  genes <- c("G1","G2","G3")
  A <- Matrix::sparseMatrix(
    i = c(1,2,2,3), j = c(2,1,3,2), x = 1,
    dims = c(3,3), dimnames = list(genes, genes)
  )
  rs <- Matrix::rowSums(A); Dinv <- Matrix::Diagonal(x = ifelse(rs > 0, 1/rs, 0))
  N <- Dinv %*% A

  mut <- matrix(0, nrow = 3, ncol = 2, dimnames = list(genes, c("P1","P2")))
  mut["G1","P1"] <- 1
  mut["G3","P2"] <- 1

  S <- rwr_multi(N, mut, alpha = 0.7, verbose = FALSE)

  testthat::expect_true(is.matrix(S))
  testthat::expect_equal(dim(S), c(length(genes), 2L))
  testthat::expect_equal(rownames(S), genes)
  testthat::expect_equal(colnames(S), c("P1","P2"))
})

# Confirms the result equals the seed vector when restart probability is 1.
testthat::test_that("alpha = 1 returns the seed vector unchanged", {
  testthat::skip_if_not_installed("Matrix")

  genes <- c("A","B","C")
  A <- Matrix::sparseMatrix(i = c(1,2,2,3), j = c(2,1,3,2), x = 1,
                            dims = c(3,3), dimnames = list(genes, genes))
  rs <- Matrix::rowSums(A); Dinv <- Matrix::Diagonal(x = ifelse(rs > 0, 1/rs, 0))
  N <- Dinv %*% A

  mut <- matrix(0, nrow = 3, ncol = 1, dimnames = list(genes, "P"))
  mut["B","P"] <- 1

  S <- rwr_multi(N, mut, alpha = 1, verbose = FALSE)

  testthat::expect_equal(drop(S[, "P"]), drop(mut[, "P"]))
})

# Checks numeric equality to the linear-system solution
testthat::test_that("solution matches closed-form linear system on a tiny graph", {
  testthat::skip_if_not_installed("Matrix")

  # G1—G2—G3, seed at G1
  genes <- c("G1","G2","G3")
  A <- Matrix::sparseMatrix(i = c(1,2,2,3), j = c(2,1,3,2), x = 1,
                            dims = c(3,3), dimnames = list(genes, genes))
  rs <- Matrix::rowSums(A); Dinv <- Matrix::Diagonal(x = ifelse(rs > 0, 1/rs, 0))
  N <- Dinv %*% A

  mut <- matrix(0, nrow = 3, ncol = 1, dimnames = list(genes, "P"))
  mut["G1","P"] <- 1

  alpha <- 0.7
  S <- rwr_multi(N, mut, alpha = alpha, eps = 1e-16, verbose = FALSE)

  # Closed-form: p = alpha * solve(I - (1 - alpha) * N, p0)
  I <- Matrix::Diagonal(n = 3)
  p0 <- mut[, "P"]
  p_expected <- as.vector(alpha * solve(as.matrix(I - (1 - alpha) * N), p0))

  # Compare numerically, ignore names
  testthat::expect_equal(unname(drop(S[, "P"])), unname(p_expected), tolerance = 1e-10)
})

# Ensures mutation rows absent from the network are treated as zeros and alignment by gene names works.
testthat::test_that("rows in mutmat not in network are treated as zeros; missing rows align", {
  testthat::skip_if_not_installed("Matrix")

  # Network genes G1,G2,G3; mutmat has G2, G4, G1 (G3 missing, G4 extra)
  genes <- c("G1","G2","G3")
  A <- Matrix::sparseMatrix(i = c(1,2,2,3), j = c(2,1,3,2), x = 1,
                            dims = c(3,3), dimnames = list(genes, genes))
  rs <- Matrix::rowSums(A); Dinv <- Matrix::Diagonal(x = ifelse(rs > 0, 1/rs, 0))
  N <- Dinv %*% A

  mut <- matrix(0, nrow = 3, ncol = 1, dimnames = list(c("G2","G4","G1"), "P"))
  mut["G4","P"] <- 1  # seed only on gene NOT in network → should be treated as zero

  S <- rwr_multi(N, mut, alpha = 0.7, verbose = FALSE)

  testthat::expect_true(all(abs(S[, "P"]) == 0))
})

# Asserts two patients with identical seeds produce identical propagated scores.
testthat::test_that("identical patients yield identical propagated profiles", {
  testthat::skip_if_not_installed("Matrix")

  genes <- c("g1","g2","g3","g4")
  A <- Matrix::sparseMatrix(
    i = c(1,2,3,4, 2,3,4,1), j = c(2,3,4,1, 1,2,3,4), x = 1,
    dims = c(4,4), dimnames = list(genes, genes)
  )
  rs <- Matrix::rowSums(A); Dinv <- Matrix::Diagonal(x = ifelse(rs > 0, 1/rs, 0))
  N <- Dinv %*% A

  mut <- matrix(0, nrow = 4, ncol = 2, dimnames = list(genes, c("P1","P2")))
  mut["g2","P1"] <- 1
  mut["g2","P2"] <- 1  # identical seed vectors

  S <- rwr_multi(N, mut, alpha = 0.8, verbose = FALSE)

  testthat::expect_equal(unname(drop(S[, "P1"])), unname(drop(S[, "P2"])), tolerance = 1e-12)
})

# Confirms a zero seed column yields an all-zero output column.
testthat::test_that("all-zero seed columns produce all-zero outputs", {
  testthat::skip_if_not_installed("Matrix")

  genes <- c("A","B")
  A <- Matrix::sparseMatrix(i = c(1,2), j = c(2,1), x = 1,
                            dims = c(2,2), dimnames = list(genes, genes))
  rs <- Matrix::rowSums(A); Dinv <- Matrix::Diagonal(x = ifelse(rs > 0, 1/rs, 0))
  N <- Dinv %*% A

  mut <- matrix(0, nrow = 2, ncol = 1, dimnames = list(genes, "P"))

  S <- rwr_multi(N, mut, alpha = 0.6, verbose = FALSE)
  testthat::expect_true(all(S == 0))
})

# Runs a minimal end-to-end computation on a small mutmat subset and tiny network to ensure nontrivial output.
testthat::test_that("rwr_multi runs on a packaged mutmat subset with a tiny network (smoke test)", {
  testthat::skip_if_not_installed("Matrix")

  # Load packaged data if available; skip otherwise
  suppressWarnings(try(data("mutmat"), silent = TRUE))
  testthat::skip_if_not(exists("mutmat"), message = "packaged data 'mutmat' not available")

  gsel <- rownames(mutmat)[1:4]
  psel <- colnames(mutmat)[1:2]

  mm <- mutmat[gsel, psel, drop = FALSE]
  mm[] <- 0
  mm[gsel[1], psel[1]] <- 1
  mm[gsel[4], psel[2]] <- 1

  # Chain 1-2-3-4 (undirected), then row-normalize
  A <- Matrix::sparseMatrix(
    i = c(1,2,3), j = c(2,3,4), x = 1,
    dims = c(4,4), dimnames = list(gsel, gsel)
  )
  A <- A + Matrix::t(A)
  rs <- Matrix::rowSums(A); Dinv <- Matrix::Diagonal(x = ifelse(rs > 0, 1/rs, 0))
  N <- Dinv %*% A

  S <- rwr_multi(N, mm, alpha = 0.7, verbose = FALSE)

  testthat::expect_equal(dim(S), dim(mm))
  testthat::expect_true(all(S >= 0))
  testthat::expect_gt(sum(S), 0)
})


testthat::test_that("rwr_multi runs on 1000 genes and 10 patients from mutmat (smoke test)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("Matrix")

  # Load packaged data; skip if not available
  suppressWarnings(try(data("mutmat"), silent = TRUE))
  testthat::skip_if_not(exists("mutmat"), message = "packaged data 'mutmat' not available")

  # Select 1000 genes and 10 patients (small subset of patients)
  gsel <- rownames(mutmat)[seq_len(min(1000L, nrow(mutmat)))]
  psel <- colnames(mutmat)[seq_len(min(10L,   ncol(mutmat)))]

  mm <- mutmat[gsel, psel, drop = FALSE]

  # Ensure non-trivial seeds: place one seed per patient at evenly spaced genes
  mm[,] <- 0
  pos <- unique(pmax(1L, round(seq(1, length(gsel), length.out = ncol(mm)))))
  for (j in seq_along(psel)) {
    mm[pos[j], j] <- 1
  }

  # Build a sparse ring network over the 1000 genes and row-normalize (row-stochastic)
  n <- length(gsel)
  i <- c(1:(n-1), 2:n, n, 1)
  j <- c(2:n,     1:(n-1), 1, n)
  A <- Matrix::sparseMatrix(i = i, j = j, x = 1,
                            dims = c(n, n), dimnames = list(gsel, gsel))
  rs <- Matrix::rowSums(A)
  N  <- Matrix::Diagonal(x = 1/rs) %*% A  # all degrees are 2 in a ring

  # Run with moderate tolerance/iterations to keep the test fast
  S <- rwr_multi(N, mm, alpha = 0.7, eps = 1e-8, maxiters = 2000, verbose = FALSE)

  # Basic sanity checks
  testthat::expect_equal(dim(S), dim(mm))
  testthat::expect_equal(rownames(S), gsel)
  testthat::expect_equal(colnames(S), psel)
  testthat::expect_true(all(is.finite(S)))
  testthat::expect_true(all(S >= 0))
  testthat::expect_gt(sum(S), 0)
  # Each patient should yield a non-zero profile
  testthat::expect_true(all(colSums(S) > 0))
})

# Topology bias check: with an all-ones seed on a row-stochastic network (no isolates),
# the steady-state should be all ones (up to rounding).
testthat::test_that("rwr_multi: no topology bias with all-ones seed on row-stochastic network", {
  testthat::skip_if_not_installed("Matrix")

  genes <- paste0("G", 1:5)
  i <- c(1:4, 2:5, 5, 1)
  j <- c(2:5, 1:4, 1, 5)
  A <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(5, 5),
                            dimnames = list(genes, genes))

  rs <- Matrix::rowSums(A)
  N  <- Matrix::Diagonal(x = 1/rs) %*% A

  mut <- matrix(1, nrow = 5, ncol = 1, dimnames = list(genes, "P"))
  S <- rwr_multi(N, mut, alpha = 0.7, eps = 1e-14, verbose = FALSE)

  # ignore names when comparing
  testthat::expect_equal(unname(round(drop(S[, "P"]), 3)), rep(1, 5))
})
