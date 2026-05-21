% Auto-generated --- do not edit by hand
\name{ctbn_fctbn}
\alias{ctbn_fctbn}
\title{Fit a Functional CTBN (Faruqui 2021)}
\usage{
ctbn_fctbn(DT_wide,
           fixed_covs = character(0),
           time_varying_covs = character(0),
           target_conditions = NULL,
           lambda = 1e-3,
           max_iter = 30000L,
           tol = 1e-7,
           gmm_stop = TRUE,
           gmm_warmup = 3000L,
           adaptive = TRUE,
           pilot_iter = 500L,
           parallel = FALSE,
           verbose = TRUE)
}
\arguments{
\item{DT_wide}{A wide-format data.table; same layout as
\code{\link{ctbn_map}}.}
\item{fixed_covs}{Character vectors of covariate
column names.}
\item{time_varying_covs}{Character vectors of covariate
column names.}
\item{target_conditions}{Optional subset of conditions to fit.}
\item{lambda}{Numeric scalar; the (non-adaptive) base regularisation
strength used to construct adaptive weights as
\eqn{\lambda_j = \lambda / \lVert \tilde\beta_j\rVert} after a
short un-penalised pilot.}
\item{max_iter}{Integer; FISTA iteration cap. Default 30 000.}
\item{tol}{Numeric; relative-change convergence tolerance.}
\item{gmm_stop}{Logical; if TRUE (default), fit a 2-component GMM to
the estimated block norms after a warm-up and zero out groups
inside the near-zero component.}
\item{gmm_warmup}{Integer; FISTA iterations before the GMM check.}
\item{adaptive}{Logical; if TRUE (default), apply adaptive group
weights from a pilot unpenalised fit.}
\item{pilot_iter}{Integer; iterations for the unpenalised pilot.}
\item{parallel}{As in \code{\link{ctbn_map}}.}
\item{verbose}{As in \code{\link{ctbn_map}}.}
}
\description{
Models the conditional intensity \eqn{q_{m|u}} as a Poisson
regression on exogenous covariates, with parent effects assumed to
be multiplicative (i.e.\ additive on the log-rate scale). Structure
learning is achieved via adaptive group lasso applied per (parent,
target) coefficient block, solved with FISTA. An optional Gaussian
mixture model early-stop snaps near-zero blocks to exactly zero,
avoiding the long convergence tail of pure proximal gradient.
}
\value{
An object of class \code{c("ctbn_fctbn_fit", "ctbn_fit")}
with the standard slots \code{beta_matrix}, \code{pip_matrix}
(1 - block-soft-threshold indicator), \code{intensity_matrix},
\code{fctbn_fits} (per-target detail) and \code{call_args}.
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
DTw <- prepare_wide(DT)

# Functional CTBN: Poisson regression with adaptive group lasso and
# multiplicative parent effects.
fit <- ctbn_fctbn(DTw, fixed_covs = c("age", "sex_male"),
                   max_iter = 5000)
summary(fit)
}

}
