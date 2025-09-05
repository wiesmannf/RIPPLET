testthat::test_that("ripplet_path equals get_npps(get_pps()) on tiny synthetic example", {
  # gene scores
  gs <- matrix(c(0.1, 0.3), nrow=2, dimnames=list(c("G1","G2"), "P1"))

  # two normalized pathway matrices (become M=-I after diag-1L)
  P1 <- matrix(0, 2, 2, dimnames=list(c("G1","G2"), c("G1","G2")))
  P2 <- matrix(0, 2, 2, dimnames=list(c("G1","G2"), c("G1","G2")))
  pw <- list(
    INFO = data.frame(),
    SPIA_Norm_Adj_Matrices = list(A = P1, B = P2),
    MaxHit = c(A = 2, B = 2)
  )

  npps1 <- ripplet_path(gene_scores = gs, pathways = pw, show_progress = FALSE)
  pps   <- get_pps(gs, pathways = pw, show_progress = FALSE)
  npps2 <- get_npps(pps, pathways = pw, show_progress = FALSE)

  testthat::expect_equal(npps1, npps2)
  testthat::expect_identical(rownames(npps1), c("A","B"))
  testthat::expect_identical(colnames(npps1), "P1")
})
