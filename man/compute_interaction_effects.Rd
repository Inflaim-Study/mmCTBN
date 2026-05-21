% Auto-generated --- do not edit by hand
\name{compute_interaction_effects}
\alias{compute_interaction_effects}
\title{Compute pairwise interaction effects (synergistic risk decomposition)}
\usage{
compute_interaction_effects(fit, DT_wide, tau = 5, pip_thresh = 0.50)
}
\arguments{
\item{fit}{A \code{ctbn_fit} object (MAP or Stan).}
\item{DT_wide}{Wide-format data.table used for fitting (needed for
covariate reference values).}
\item{tau}{Numeric, time horizon (years).}
\item{pip_thresh}{Numeric in [0, 1]; significance threshold for
interaction selection.}
}
\description{
For each significant triplet \eqn{(j, k \to m)} where the interaction
PIP exceeds \code{pip_thresh}, returns the parameter decomposition
(main RRs, synergistic multiplier, joint RR) and the absolute
synergistic excess
\eqn{\Delta F_m(\tau) = F(j,k) - F(j) - F(k) + F(\emptyset)} at
time horizon \code{tau}.
}
\value{
A data.table with columns:
\describe{
\item{target_condition, condition_j, condition_k}{names}
\item{beta_j, beta_k, beta_jk}{coefficients}
\item{RR_j, RR_k, syn_multiplier, joint_RR}{rate ratios}
\item{F_baseline, F_j_only, F_k_only, F_joint}{absolute risks}
\item{delta_F}{synergistic excess}
\item{pip_jk}{interaction selection score}
}
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500,
                           interaction_order = 2, seed = 1)
DTw <- prepare_wide(DT)

# Interaction decomposition needs a model with max_order >= 2.
fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 2)
syn <- compute_interaction_effects(fit, DTw, tau = 5, pip_thresh = 0.5)
head(syn)
}

}
