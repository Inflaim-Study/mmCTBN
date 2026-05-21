skip_if_no_methods_data <- function() {
  if (!exists(".methods_DT", envir = globalenv())) {
    nw <- make_random_network(n_nodes = 5, seed = 11)
    DT <- simulate_ctbn_data(network = nw,
                              n_patients = 120,
                              t_horizon = 6,
                              interaction_order = 2,
                              seed = 11)
    assign(".methods_DT", DT, envir = globalenv())
    assign(".methods_NW", nw, envir = globalenv())
  }
  get(".methods_DT", envir = globalenv())
}

test_that("ctbn_classic produces a ctbn_classic_fit with expected slots", {
  DT <- skip_if_no_methods_data()
  fit <- ctbn_classic(DT,
                       fixed_covs        = c("age", "sex_male"),
                       time_varying_covs = c("smk_current", "smk_former"),
                       max_parents       = 2L,
                       parallel          = FALSE,
                       verbose           = FALSE)
  expect_s3_class(fit, "ctbn_classic_fit")
  expect_s3_class(fit, "ctbn_fit")
  expect_true(is.matrix(fit$beta_matrix))
  expect_equal(rownames(fit$beta_matrix), colnames(fit$beta_matrix))
  expect_true("rates" %in% names(fit))
  # Each target should have a 2-column rates matrix
  for (tgt in names(fit$rates)) {
    rm <- fit$rates[[tgt]]
    expect_true(is.matrix(rm))
    expect_equal(colnames(rm), c("lambda_01", "lambda_10"))
    expect_true(all(rm >= 0, na.rm = TRUE))
  }
})

test_that("get_lp.ctbn_classic_fit returns sensible lp/lambda", {
  DT  <- skip_if_no_methods_data()
  fit <- ctbn_classic(DT, max_parents = 2L,
                       parallel = FALSE, verbose = FALSE)
  tgt <- colnames(fit$beta_matrix)[1]
  res <- get_lp(fit, DT, tgt)
  expect_named(res, c("lp", "lambda"))
  expect_equal(length(res$lp), nrow(DT))
  # Patients already in the target state should return NA
  is_one <- as.integer(DT[[tgt]] == 1L)
  expect_true(all(is.na(res$lambda[is_one == 1L])))
  # Patients in state 0 should have non-NA, non-negative lambdas
  expect_true(all(is.finite(res$lambda[is_one == 0L])))
  expect_true(all(res$lambda[is_one == 0L] >= 0))
})

test_that("ctbn_fctbn produces a ctbn_fctbn_fit and selects a sensible edge set", {
  DT <- skip_if_no_methods_data()
  fit <- ctbn_fctbn(DT,
                     fixed_covs        = c("age", "sex_male"),
                     time_varying_covs = c("smk_current", "smk_former"),
                     lambda            = 0.05,
                     max_iter          = 1500L,
                     gmm_warmup        = 500L,
                     pilot_iter        = 200L,
                     adaptive          = TRUE,
                     parallel          = FALSE,
                     verbose           = FALSE)
  expect_s3_class(fit, "ctbn_fctbn_fit")
  expect_s3_class(fit, "ctbn_fit")
  expect_true(is.matrix(fit$beta_matrix))
  expect_true(all(is.finite(fit$beta_matrix)))
  # pip_matrix entries are 0/1/NA
  pm <- fit$pip_matrix
  expect_true(all(pm %in% c(0, 1) | is.na(pm)))
})

test_that("get_lp.ctbn_fctbn_fit returns finite predictions", {
  DT  <- skip_if_no_methods_data()
  fit <- ctbn_fctbn(DT,
                     fixed_covs = c("age"),
                     lambda     = 0.05,
                     max_iter   = 1000L,
                     pilot_iter = 100L,
                     parallel   = FALSE, verbose = FALSE)
  tgt <- colnames(fit$beta_matrix)[1]
  res <- get_lp(fit, DT, tgt)
  expect_named(res, c("lp", "lambda"))
  expect_equal(length(res$lp), nrow(DT))
  expect_true(all(is.finite(res$lambda)))
})

test_that("ctbn_cph returns a ctbn_cph_fit and runs end-to-end", {
  skip_if_not_installed("survival")
  DT <- skip_if_no_methods_data()
  fit <- ctbn_cph(DT,
                   fixed_covs        = c("age", "sex_male"),
                   time_varying_covs = c("smk_current"),
                   max_parents       = 2L,
                   alpha             = 0.20,        # permissive for a small sim
                   parallel          = FALSE,
                   verbose           = FALSE)
  expect_s3_class(fit, "ctbn_cph_fit")
  expect_s3_class(fit, "ctbn_fit")
  expect_true("cox_fits" %in% names(fit))
  expect_true("beta_cov" %in% names(fit))
  expect_s3_class(fit$classic, "ctbn_classic_fit")
})

test_that("get_lp.ctbn_cph_fit produces individualised lambdas", {
  skip_if_not_installed("survival")
  DT  <- skip_if_no_methods_data()
  fit <- ctbn_cph(DT,
                   fixed_covs  = c("age", "sex_male"),
                   max_parents = 2L,
                   alpha       = 0.20,
                   parallel    = FALSE,
                   verbose     = FALSE)
  tgt <- colnames(fit$beta_matrix)[1]
  res <- get_lp(fit, DT, tgt)
  expect_equal(length(res$lp), nrow(DT))
  # At least some rows must be active (target = 0); they should be finite
  is_zero <- DT[[tgt]] == 0L
  expect_true(any(is_zero))
  expect_true(all(is.finite(res$lambda[is_zero])))
})

test_that("ctbn_fit dispatcher reaches every method", {
  DT <- skip_if_no_methods_data()

  method_args <- list(
    map     = list(prior = "spike_slab", max_order = 1L, compute_se = FALSE),
    classic = list(max_parents = 2L),
    fctbn   = list(max_iter = 800L, pilot_iter = 100L, lambda = 0.05))

  for (m in c("map", "classic", "fctbn")) {
    fit <- do.call(ctbn_fit, c(
      list(DT_wide = DT, method = m,
           fixed_covs = c("age"),
           verbose = FALSE),
      method_args[[m]]))
    expect_s3_class(fit, "ctbn_fit")
  }

  if (requireNamespace("survival", quietly = TRUE)) {
    fit_cph <- ctbn_fit(DT, method = "cph",
                         fixed_covs  = c("age"),
                         max_parents = 2L,
                         alpha       = 0.20,
                         verbose     = FALSE)
    expect_s3_class(fit_cph, "ctbn_cph_fit")
  }
})
