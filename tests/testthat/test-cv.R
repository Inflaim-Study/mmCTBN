test_that("make_patient_folds returns disjoint patient assignments", {
  DT  <- data.table::data.table(eid = rep(1:20, each = 3),
                                 time_to_event = rep(1:3, 20),
                                 dt  = 1)
  fld <- make_patient_folds(DT, k = 4, seed = 0)
  expect_equal(nrow(fld), 20)
  expect_setequal(fld$fold, 1:4)
})

test_that("summarise_cv handles an empty AUC slice gracefully", {
  cv_dt <- data.table::data.table(
    fold        = rep(1:2, each = 4),
    model       = rep(c("a", "b"), times = 4),
    target      = rep(c("X", "Y"), each = 2, times = 2),
    poisson_ll  = rnorm(8),
    brier       = runif(8),
    eval_time   = NA_real_,
    tdauc       = NA_real_,
    n_test_rows = 100L,
    n_test_pats = 50L)
  s <- summarise_cv(cv_dt)
  expect_named(s, c("scalar", "tdauc"))
  expect_equal(nrow(s$tdauc), 0)
  expect_true(all(c("mean_poisson_ll", "mean_brier") %in% names(s$scalar)))
})
