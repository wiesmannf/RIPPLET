testthat::test_that("impact_scoring: alignment, extras dropped, finite [0,1], ranking preserved", {
  R <- matrix(c(1,0,0, 1,0,0, 0,1,0), nrow=3,
              dimnames=list(paste0("F",1:3), c("P1","P2","P3")))
  S <- diag(1,3); rownames(S) <- colnames(S) <- c("P1","P2","P3")
  S["P1","P2"] <- S["P2","P1"] <- 1  # strong mutual similarity
  out <- impact_scoring(R[,1:2,drop=FALSE], S, delta = 0.2)
  testthat::expect_equal(dim(out), c(3L,2L))
  testthat::expect_identical(colnames(out), c("P1","P2"))
  testthat::expect_true(all(is.finite(out)))
  testthat::expect_true(all(out >= 0 & out <= 1))
  # F1 is maximal for both since P1,P2 share the same profile on F1
  testthat::expect_equal(rownames(out)[which.max(out[, "P1"])], "F1")
  testthat::expect_equal(rownames(out)[which.max(out[, "P2"])], "F1")
})

testthat::test_that("impact_scoring: error when sim_mat missing a patient", {
  R <- matrix(c(1,0,0, 0,1,0), nrow=3,
              dimnames=list(paste0("F",1:3), c("P1","P2")))
  S <- diag(1,1); rownames(S) <- colnames(S) <- "P1"
  testthat::expect_error(impact_scoring(R, S), "missing similarity data")
})

testthat::test_that("impact_scoring: extreme delta values still valid", {
  R <- matrix(runif(12), nrow=4, dimnames=list(paste0("F",1:4), c("A","B","C")))
  S <- matrix(1, 3, 3); diag(S) <- 1; rownames(S) <- colnames(S) <- c("A","B","C")
  o1 <- impact_scoring(R, S, delta = 1e-6)
  o2 <- impact_scoring(R, S, delta = 10)
  testthat::expect_true(all(is.finite(o1)) && all(is.finite(o2)))
  testthat::expect_true(all(o1 >= 0 & o1 <= 1) && all(o2 >= 0 & o2 <= 1))
})
