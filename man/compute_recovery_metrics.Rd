% Auto-generated --- do not edit by hand
\name{compute_recovery_metrics}
\alias{compute_recovery_metrics}
\title{Recovery metrics (bias, RMSE, coverage of credible intervals)}
\usage{
compute_recovery_metrics(fit, network)
}
\arguments{
\item{fit}{A fitted \code{ctbn_fit} (must have SEs -- typically a MAP
fit with \code{compute_se = TRUE}).}
\item{network}{The true network spec.}
}
\description{
Recovery metrics (bias, RMSE, coverage of credible intervals)
}
\value{
A data.table with rows per parameter and columns
\code{condition}, \code{parameter}, \code{param_type}, \code{true_val},
\code{est_mean}, \code{ci_lo}, \code{ci_hi}, \code{bias},
\code{sq_err}, \code{covered}.
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
DTw <- prepare_wide(DT)

# Recovery metrics need standard errors (compute_se = TRUE).
fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 1,
                 compute_se = TRUE)
rec <- compute_recovery_metrics(fit, network = net)
head(rec)
}

}
