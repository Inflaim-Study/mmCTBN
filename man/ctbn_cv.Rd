% Auto-generated --- do not edit by hand
\name{ctbn_cv}
\alias{ctbn_cv}
\title{Patient-level k-fold cross-validation for CTBN models}
\usage{
ctbn_cv(DT_wide,
        fit_fns,
        fixed_covs = character(0),
        time_varying_covs = character(0),
        k_folds = 3,
        targets = NULL,
        seed = 42,
        eval_times = NULL,
        .progress = FALSE,
        ...)
}
\arguments{
\item{DT_wide}{Wide-format data.table (one row per at-risk interval
per patient). Must contain \code{eid}, \code{time_to_event},
\code{dt}, the condition columns and any covariates.}
\item{fit_fns}{Named list of fitting functions. Each must take a
data.table as its first argument and accept
\code{target_conditions} as a named argument. Most commonly these
are wrappers around \code{ctbn_map} and/or \code{ctbn_stan}.}
\item{fixed_covs}{Character vector of fixed covariate names.}
\item{time_varying_covs}{Character vector of time-varying covariate
names.}
\item{k_folds}{Integer, number of folds (default 3).}
\item{targets}{Optional character subset of conditions to evaluate.}
\item{seed}{Integer RNG seed for fold assignment and worker streams.}
\item{eval_times}{Optional numeric vector of evaluation times for
time-dependent AUC (requires \pkg{timeROC}).}
\item{.progress}{Logical; if TRUE and \pkg{progressr} is installed,
show progress.}
\item{...}{Additional arguments forwarded to every fit function.}
}
\description{
Compares any number of CTBN model specifications under a common k-fold
patient-level split. Per-fold computation is dispatched via
\code{future.apply::future_lapply}, so set the desired plan() before
calling, e.g.\ \code{future::plan(multisession, workers = 4)}.
}
\value{
A long-format data.table with columns: \code{fold, model,
target, poisson_ll, brier, eval_time, tdauc, n_test_rows,
n_test_pats}.
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 400,
                           interaction_order = 2, seed = 1)
DTw <- prepare_wide(DT)

# Compare a MAP fit against a classical CTBN under 3-fold CV.
# Each fit function takes the data as its first argument and must
# accept `target_conditions`; ctbn_cv() injects it per fold.
fit_fns <- list(
  map     = function(DT, ...) ctbn_map(DT, prior = "spike_slab",
                                         max_order = 2, ...),
  classic = function(DT, ...) ctbn_classic(DT, max_parents = 3, ...))

future::plan(future::sequential)
cv <- ctbn_cv(DTw, fit_fns = fit_fns, k_folds = 3,
               eval_times = c(2, 5), seed = 1)
head(cv)
summarise_cv(cv)
}

}
