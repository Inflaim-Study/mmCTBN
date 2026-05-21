test_that("make_random_network runs for a range of node counts", {
  for (p in c(3, 5, 8)) {
    nw <- make_random_network(n_nodes = p, seed = 1,
                                max_parents = min(4L, p - 1L))
    expect_equal(length(nw$conds), p)
    expect_equal(length(nw$beta0), p)
    expect_true(is.list(nw$beta_int2))
  }
})

test_that("make_default_network returns a 10-node network", {
  nw <- make_default_network()
  expect_equal(length(nw$conds), 10)
  expect_true(all(c("HTN", "DM", "CLRD") %in% nw$conds))
})

test_that("simulate_ctbn_data respects interaction_order in [1, 5]", {
  for (k in 1:3) {
    nw <- make_random_network(n_nodes = 5, seed = 42)
    DT <- simulate_ctbn_data(network = nw,
                              n_patients = 30,
                              t_horizon = 5,
                              interaction_order = k,
                              seed = 7)
    expect_s3_class(DT, "data.table")
    expect_true(all(c("eid", "time_to_event", "dt") %in% names(DT)))
    expect_true(nrow(DT) > 0)
  }
})

test_that("simulate_ctbn_data rejects interaction_order outside [1,5]", {
  nw <- make_random_network(n_nodes = 3, seed = 1, max_parents = 2L)
  expect_error(simulate_ctbn_data(network = nw,
                                   n_patients = 5,
                                   interaction_order = 0L),
                "interaction_order")
  expect_error(simulate_ctbn_data(network = nw,
                                   n_patients = 5,
                                   interaction_order = 6L),
                "interaction_order")
})
