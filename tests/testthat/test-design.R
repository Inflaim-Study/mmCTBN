test_that("build_interaction_cols returns expected number of columns", {
  d <- data.frame(A = c(0, 1, 0, 1),
                  B = c(0, 0, 1, 1),
                  C = c(1, 1, 0, 0))
  res <- build_interaction_cols(d, c("A", "B", "C"), k = 2)
  # All pairwise products: AB, AC, BC = 3
  expect_equal(length(res$cols), 3)
  # 2-way interactions are encoded as order = 1 (ord - 1L)
  expect_equal(unique(res$orders), 1L)

  res3 <- build_interaction_cols(d, c("A", "B", "C"), k = 3)
  # 3 pairs + 1 triplet = 4
  expect_equal(length(res3$cols), 4)
  expect_true(any(res3$orders == 2L))
})

test_that("build_interaction_cols caps at k = 5", {
  d <- as.data.frame(matrix(rbinom(6 * 7, 1, 0.5), 6, 7,
                             dimnames = list(NULL, LETTERS[1:7])))
  expect_error(build_interaction_cols(d, LETTERS[1:7], k = 6), "k")
})
