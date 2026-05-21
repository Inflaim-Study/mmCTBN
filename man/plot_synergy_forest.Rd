% Auto-generated --- do not edit by hand
\name{plot_synergy_forest}
\alias{plot_synergy_forest}
\title{Forest plot of synergistic interaction effects}
\usage{
plot_synergy_forest(synergy_dt, top_n = 25, x = c("beta_jk", "joint_log_RR"))
}
\arguments{
\item{synergy_dt}{A data.table from \code{compute_interaction_effects}.}
\item{top_n}{Integer, number of effects to display.}
\item{x}{One of \code{"beta_jk"} (interaction coef, default) or
\code{"joint_log_RR"} (sum of mains + interaction).}
}
\description{
Plots \code{beta_jk} (or the joint log-RR) for the top
\code{top_n} interactions by absolute synergistic excess from
\code{\link{compute_interaction_effects}}.
}
\value{
A \code{ggplot} object.
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500,
                           interaction_order = 2, seed = 1)
DTw <- prepare_wide(DT)
fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 2)

syn <- compute_interaction_effects(fit, DTw, tau = 5)
plot_synergy_forest(syn, top_n = 15)
}

}
