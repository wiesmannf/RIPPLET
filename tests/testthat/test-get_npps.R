# tests/testthat/test-get_npps.R

testthat::test_that("get_npps: basic normalization with named MaxHit vector", {
  gs  <- matrix(c(0.005, 0.02), nrow = 2,
                dimnames = list(c("G1","G2"), "P1"))
  P_ok <- matrix(0, 2, 2, dimnames = list(c("G1","G2"), c("G1","G2")))  # M = -I
  pws  <- list(P_ok = P_ok)

  pps <- get_pps(gene_scores = gs, pathways = pws,
                 threshold_level_ps = 0, show_progress = FALSE)

  mh_vec <- c(P_ok = sum(gs[, "P1"]))  # so normalized score is exactly 1

  npps <- get_npps(pps = pps, pathways = mh_vec, show_progress = FALSE)

  testthat::expect_true(is.matrix(npps))
  testthat::expect_identical(rownames(npps), names(mh_vec))
  testthat::expect_identical(colnames(npps), "P1")
  testthat::expect_equal(as.numeric(npps["P_ok", "P1"]), 1)
})

testthat::test_that("get_npps: accepts full pathways list with MaxHit; drops non-overlaps", {
  gs  <- matrix(c(0.1, 0.2), nrow = 2,
                dimnames = list(c("G1","G2"), "P1"))
  P_ok <- matrix(0, 2, 2, dimnames = list(c("G1","G2"), c("G1","G2")))
  full_pw <- list(
    SPIA_Norm_Adj_Matrices = list(P_ok = P_ok),
    MaxHit = c(P_ok = sum(gs[, "P1"]), P_extra = 123) # extra will be dropped
  )

  pps <- get_pps(gene_scores = gs, pathways = full_pw,
                 threshold_level_ps = 0, show_progress = FALSE)
  npps <- get_npps(pps = pps, pathways = full_pw, show_progress = FALSE)

  testthat::expect_true("P_ok" %in% rownames(npps))
  testthat::expect_equal(as.numeric(npps["P_ok", "P1"]), 1)
  testthat::expect_false("P_extra" %in% rownames(npps))
})

testthat::test_that("get_npps: pps validation errors (empty, unnamed, missing Pathway)", {
  # empty pps
  testthat::expect_error(get_npps(list(), pathways = c(A = 1)), "non-empty")

  # pps must be a named list
  unnamed <- list(data.frame(Pathway = "A", sumpfs = 1))
  testthat::expect_error(get_npps(unnamed, pathways = c(A = 1)), "named list")

  # first element missing Pathway column triggers error
  bad1 <- list(P1 = data.frame(NotPathway = "A", sumpfs = 1))
  testthat::expect_error(get_npps(bad1, pathways = c(A = 1)), "Pathway")
})

testthat::test_that("get_npps: missing 'sumpfs' column yields NAs (no error)", {
  # Function currently treats missing sumpfs softly (fills NA), not as an error
  good_structure <- list(P1 = data.frame(Pathway = "A"))
  npps <- get_npps(good_structure, pathways = c(A = 2), show_progress = FALSE)

  testthat::expect_true(is.matrix(npps))
  testthat::expect_identical(rownames(npps), "A")
  testthat::expect_identical(colnames(npps), "P1")
  testthat::expect_true(is.na(npps["A", "P1"]))
})

testthat::test_that("get_npps: pathways argument validation and overlap checks", {
  good <- list(P1 = data.frame(Pathway = "A", sumpfs = 1))

  # wrong type (numeric without names) → generic type error
  testthat::expect_error(get_npps(good, pathways = 1), regexp = NULL)

  # empty MaxHit in full list → explicit 'empty' error
  testthat::expect_error(get_npps(good, pathways = list(MaxHit = numeric(0))), "empty")

  # no overlap between pps pathways and MaxHit names → explicit overlap error
  testthat::expect_error(get_npps(good, pathways = c(B = 1)), "No overlapping pathways")
})

testthat::test_that("get_npps: zero/negative MaxHit produce NAs", {
  pps <- list(
    P1 = data.frame(Pathway = c("A","B"), sumpfs = c(2, 3))
  )
  mh  <- c(A = 0, B = -5)  # both invalid denominators

  npps <- get_npps(pps, pathways = mh, show_progress = FALSE)

  testthat::expect_true(all(is.na(npps[, "P1"])))
})

testthat::test_that("get_npps: multiple samples align correctly and subset to overlap", {
  pps <- list(
    P1 = data.frame(Pathway = c("A","B"), sumpfs = c(2, 3)),
    P2 = data.frame(Pathway = c("B","C"), sumpfs = c(4, 5))
  )
  mh <- c(A = 2, B = 4)  # C is absent → will be dropped

  npps <- get_npps(pps, pathways = mh, show_progress = FALSE)

  testthat::expect_identical(rownames(npps), c("A","B"))
  testthat::expect_identical(colnames(npps), c("P1","P2"))
  # P1: A=2/2=1, B=3/4=0.75
  testthat::expect_equal(npps["A","P1"], 1)
  testthat::expect_equal(npps["B","P1"], 0.75)
  # P2: A missing -> NA, B=4/4=1
  testthat::expect_true(is.na(npps["A","P2"]))
  testthat::expect_equal(npps["B","P2"], 1)
})
