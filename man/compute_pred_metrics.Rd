% Auto-generated --- do not edit by hand
\name{compute_pred_metrics}
\alias{compute_pred_metrics}
\title{Predictive metrics on test data (Poisson log-likelihood, Brier, AUC)}
\usage{
compute_pred_metrics(fit, DT_test, eval_times = 1L)
}
\arguments{
\item{fit}{A fitted \code{ctbn_fit}.}
\item{DT_test}{Test data.}
\item{eval_times}{Numeric vector of evaluation horizons for
time-dependent AUC (default \code{1L}).}
}
\description{
Predictive metrics on test data (Poisson log-likelihood, Brier, AUC)
}
\value{
A data.table with columns \code{condition}, \code{pll},
\code{brier}, plus one column per \code{eval_times}.
}
\examples{
\dontrun{
net    <- make_random_network(n_nodes = 6, seed = 1)
DT     <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
DTw    <- prepare_wide(DT)
DTtest <- prepare_wide(
  simulate_ctbn_data(network = net, n_patients = 300, seed = 2))

fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 1)
compute_pred_metrics(fit, DTtest, eval_times = c(2, 5))
}

}
