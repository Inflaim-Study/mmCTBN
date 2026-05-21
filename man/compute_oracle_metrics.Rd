% Auto-generated --- do not edit by hand
\name{compute_oracle_metrics}
\alias{compute_oracle_metrics}
\title{Oracle (true-rate) predictive metrics on test data}
\usage{
compute_oracle_metrics(DT_test, network, interaction_order = 3L)
}
\arguments{
\item{DT_test}{A wide-format data.table including the per-target
event indicators (\code{<cond>_event}).}
\item{network}{The network spec used to generate the data.}
\item{interaction_order}{Integer in [1, 5]; truth's interaction order.}
}
\description{
Computes the Poisson log-likelihood and Brier score using the
\emph{true} simulation rates -- i.e., the upper-bound performance
achievable if the data-generating process were known.
}
\value{
A data.table with columns \code{condition}, \code{oracle_pll},
\code{oracle_brier}.
}
\examples{
\dontrun{
net    <- make_random_network(n_nodes = 6, seed = 1)
DTtest <- prepare_wide(
  simulate_ctbn_data(network = net, n_patients = 300,
                      interaction_order = 2, seed = 2))

# Upper-bound performance under the known data-generating process.
compute_oracle_metrics(DTtest, network = net, interaction_order = 2)
}

}
