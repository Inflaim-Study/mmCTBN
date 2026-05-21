% Auto-generated --- do not edit by hand
\name{summarise_cv}
\alias{summarise_cv}
\title{Summarise cross-validation results}
\usage{
summarise_cv(cv_results)
}
\arguments{
\item{cv_results}{A data.table from \code{ctbn_cv}.}
}
\description{
Aggregates fold-level metrics from \code{\link{ctbn_cv}}.
}
\value{
A list with two elements:
\describe{
\item{scalar}{data.table with mean/SE of \code{poisson_ll} and
\code{brier} per (model, target).}
\item{tdauc}{data.table with mean/SE of time-dependent AUC per
(model, target, eval_time).}
}
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 400, seed = 1)
DTw <- prepare_wide(DT)

fit_fns <- list(map = function(DT, ...) ctbn_map(DT, max_order = 1, ...))
cv      <- ctbn_cv(DTw, fit_fns = fit_fns, k_folds = 3,
                    eval_times = c(2, 5), seed = 1)

summ <- summarise_cv(cv)
summ$scalar   # mean / SE of Poisson log-lik and Brier per model+target
summ$tdauc    # mean / SE of time-dependent AUC
}

}
