% Auto-generated --- do not edit by hand
\name{ctbn_classic}
\alias{ctbn_classic}
\title{Fit a classical Continuous-Time Bayesian Network (Nodelman 2002)}
\usage{
ctbn_classic(DT_wide,
             parents = NULL,
             max_parents = 4L,
             fixed_covs = character(0),
             time_varying_covs = character(0),
             target_conditions = NULL,
             add_laplace = 0.5,
             parallel = FALSE,
             verbose = TRUE)
}
\arguments{
\item{DT_wide}{A wide-format data.table; same layout as for
\code{\link{ctbn_map}}.}
\item{parents}{Optional named list: \code{parents[[m]]} is the
character vector of parent condition names for target \code{m}. If
\code{NULL}, every other condition is treated as a parent (the
fully connected CTBN).}
\item{max_parents}{Integer; cap on \code{|parents[[m]]|}. The
parameter count per target is \eqn{2^{|parents|+1}}, so values
above 6 are usually impractical.}
\item{fixed_covs}{Currently ignored (classical
CTBNs have no covariate channel). Retained in the signature so
that \code{ctbn_cv()} drivers can pass them uniformly across
methods.}
\item{time_varying_covs}{Currently ignored (classical
CTBNs have no covariate channel). Retained in the signature so
that \code{ctbn_cv()} drivers can pass them uniformly across
methods.}
\item{target_conditions}{Optional subset of conditions to fit as
targets.}
\item{add_laplace}{Numeric pseudo-count added to both
\eqn{M[\cdot]} and \eqn{T[\cdot]} to stabilise rate estimates for
parent configurations that are rare or unobserved. Default 0.5.}
\item{parallel}{As in \code{\link{ctbn_map}}.}
\item{verbose}{As in \code{\link{ctbn_map}}.}
}
\description{
Fully parametric CTBN with discrete (0/1) state variables and no
covariates. For every target node and every joint configuration of
its parents the 0->1 and 1->0 transition rates are estimated by
closed-form maximum-likelihood:
\deqn{\hat{\lambda}_{m\,|\,u} \;=\; M[m\,|\,u] \,/\, T[m\,|\,u].}

This is the original Nodelman formulation. It does not condition the
rates on exogenous covariates and grows exponentially in the number
of parents per node — keep \code{max_parents} small.
}
\value{
An object of class \code{c("ctbn_classic_fit", "ctbn_fit")}
with elements
\describe{
\item{rates}{Named list \code{rates[[m]]} = numeric matrix with
one row per parent configuration and columns
\code{lambda_01}, \code{lambda_10}.}
\item{parent_index}{Named list mapping parent-configuration
strings (e.g.\ "0_1_0") to row indices.}
\item{parents}{Echo of the \code{parents} argument.}
\item{beta_matrix}{Influencer-by-target log-RR matrix
(constructed by marginalising parent-configuration rates so
that the rest of the package — plot, summary, metrics — keeps
working).}
\item{pip_matrix, kappa_matrix}{Inclusion proxies; set to 1 for
every (parent, target) edge in the supplied graph and NA
elsewhere.}
\item{call_args}{Echo of the call.}
}
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
DTw <- prepare_wide(DT)

# Classical CTBN: fully parametric, no covariates, parents selected
# up to max_parents per node.
fit <- ctbn_classic(DTw, max_parents = 3)
summary(fit)
coef(fit, "intensity")
}

}
