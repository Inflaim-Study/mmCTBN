skip_if_no_data <- function() {
  if (!exists(".sim_DT", envir = globalenv())) {
    nw <- make_random_network(n_nodes = 5, seed = 1)
    DT <- simulate_ctbn_data(network = nw,
                              n_patients = 80,
                              t_horizon = 6,
                              interaction_order = 2,
                              seed = 1)
    assign(".sim_DT", DT, envir = globalenv())
    assign(".sim_NW", nw, envir = globalenv())
  }
  get(".sim_DT", envir = globalenv())
}

test_that("ctbn_map fits and returns the expected structure", {
  DT <- skip_if_no_data()
  fit <- ctbn_map(DT, prior = "spike_slab",
                  max_order = 2, parallel = FALSE,
                  compute_se = FALSE,
                  verbose = FALSE)

  expect_s3_class(fit, "ctbn_map_fit")
  expect_s3_class(fit, "ctbn_fit")
  expect_true(is.matrix(fit$beta_matrix))
  expect_true(is.matrix(fit$pip_matrix))
  expect_equal(rownames(fit$beta_matrix), colnames(fit$beta_matrix))
})

test_that("get_lp returns vectors with the correct length", {
  DT  <- skip_if_no_data()
  fit <- ctbn_map(DT, prior = "lasso",
                  max_order = 2, parallel = FALSE,
                  compute_se = FALSE,
                  verbose = FALSE)

  tgt <- colnames(fit$beta_matrix)[1]
  res <- get_lp(fit, DT, tgt)
  expect_named(res, c("lp", "lambda"))
  expect_equal(length(res$lp), nrow(DT))
  expect_true(all(is.finite(res$lambda) | is.na(res$lambda)))
})

test_that("ctbn_fit dispatcher returns the same class as ctbn_map", {
  DT  <- skip_if_no_data()
  fit <- ctbn_fit(DT, method = "map", prior = "horseshoe",
                  max_order = 2, parallel = FALSE,
                  compute_se = FALSE,
                  verbose = FALSE)
  expect_s3_class(fit, "ctbn_map_fit")
})
