testthat::test_that("get_pps: input validation", {
  gs <- matrix(1, 2, 1)
  colnames(gs) <- "P1"
  testthat::expect_error(get_pps(gs, pathways = list()), "row names")
  rownames(gs) <- c("G1","G2")
  testthat::expect_error(get_pps(gs, pathways = list()), "appears to be empty")
})

testthat::test_that("get_pps: solves invertible, NA for singular; thresholding works", {
  gs <- matrix(c(0.005, 0.02), nrow = 2,
               dimnames = list(c("G1","G2"), "P1"))

  # Two pathways: P_ok (invertible after shift), P_sing (singular)
  P_ok   <- matrix(0, 2, 2, dimnames = list(c("G1","G2"), c("G1","G2")))   # M = -I
  P_sing <- diag(1, 2); dimnames(P_sing) <- list(c("G1","G2"), c("G1","G2"))  # M = 0
  pw <- list(P_ok = P_ok, P_sing = P_sing)

  res1 <- get_pps(gs, pathways = pw, threshold_level_ps = 0, show_progress = FALSE)
  testthat::expect_true(is.data.frame(res1$P1))
  testthat::expect_true("sumpfs" %in% names(res1$P1))
  testthat::expect_true(is.na(res1$P1$sumpfs[res1$P1$Pathway == "P_sing"]))

  # For P_ok, M = -I → pf = X; sumpfs = sum(gs)
  sm_ok <- res1$P1$sumpfs[res1$P1$Pathway == "P_ok"]
  testthat::expect_equal(as.numeric(sm_ok), sum(gs[, "P1"]))

  # Threshold removes G1 contribution (0.005 < 0.01)
  res2 <- get_pps(gs, pathways = pw, threshold_level_ps = 0.01, show_progress = FALSE)
  sm_ok2 <- res2$P1$sumpfs[res2$P1$Pathway == "P_ok"]
  testthat::expect_equal(as.numeric(sm_ok2), 0.02)
})

