% Auto-generated --- do not edit by hand
\name{ctbn_stan}
\alias{ctbn_stan}
\title{Fit a CTBN by full Bayesian MCMC (Stan)}
\usage{
ctbn_stan(DT_wide,
          prior = c("spike_slab", "structured", "lasso", "horseshoe"),
          max_order = 1L,
          fixed_covs = character(0),
          time_varying_covs = character(0),
          target_conditions = NULL,
          variable_select = FALSE,
          pip_threshold = 0.5,
          theta = 1.0,
          a0 = 2.0,
          b0 = 1.0,
          pi0 = 0.5,
          spike_var = 0.01,
          tau0 = 1.0,
          slab_df = 4.0,
          slab_scale = 2.0,
          chains = 4,
          iter = 2000,
          warmup = 1000,
          seed = 42L,
          parallel = FALSE,
          verbose = TRUE)
}
\arguments{
\item{DT_wide}{A data.table in wide format with columns:
\itemize{
\item \code{eid}            -- patient identifier
\item \code{time_to_event}  -- start of at-risk interval
\item \code{dt}             -- length of at-risk interval
\item one column per condition (0/1 indicator)
\item one column per fixed and time-varying covariate
}}
\item{prior}{Character; one of \code{"spike_slab"}, \code{"structured"},
\code{"lasso"}, \code{"horseshoe"}.}
\item{max_order}{Integer in [1, 5]; maximum interaction order
(1 = main effects only, 5 = up to 5-way interactions).}
\item{fixed_covs}{Character vector of fixed-covariate column names.}
\item{time_varying_covs}{Character vector of time-varying covariate names.}
\item{target_conditions}{Optional character subset of conditions to fit
as targets. Default \code{NULL} = fit all.}
\item{variable_select}{Logical; if \code{TRUE}, gate the intensity matrix
by the selection threshold \code{pip_threshold}.}
\item{pip_threshold}{Numeric in [0, 1]; default 0.5.}
\item{theta}{Numeric, order penalty strength (>= 0).}

\item{a0}{Numeric (> 0), shape hyperparameter of the inverse-gamma prior
placed on the slab variance.}

\item{b0}{Numeric (> 0), rate (scale) hyperparameter of the inverse-gamma
prior placed on the slab variance.}

\item{pi0}{Base inclusion probability (spike_slab only).}
\item{spike_var}{Spike variance (spike_slab only).}
\item{tau0}{Half-Cauchy global scale (horseshoe only).}
\item{slab_df}{Slab degrees of freedom (horseshoe only).}
\item{slab_scale}{Slab scale (horseshoe only).}
\item{chains}{Integer, number of MCMC chains (default 4).}
\item{iter}{Integer, total iterations per chain (default 2000).}
\item{warmup}{Integer, warmup iterations (default 1000).}
\item{seed}{Integer, random seed.}
\item{parallel}{Logical, run nodes in parallel via \pkg{future.apply}.}
\item{verbose}{Logical, print progress messages.}
}
\description{
Fits one Stan model per target node, using Hamiltonian Monte Carlo
(NUTS) to sample from the joint posterior over influencer effects,
interactions and covariate coefficients.
}
\value{
An object of class \code{c("ctbn_stan_fit", "ctbn_fit")} with
slots \code{beta_matrix}, \code{pip_matrix}, \code{kappa_matrix},
\code{intensity_matrix}, \code{stan_fits}, \code{call_args}, ...
}
\examples{
\dontrun{
net <- make_default_network()
DT  <- simulate_ctbn_data(network = net, n_patients = 1000,
interaction_order = 2, seed = 1)
DTw <- prepare_wide(DT)

fit <- ctbn_stan(DTw, prior = "horseshoe", max_order = 2,
fixed_covs = c("intercept","sex_male"),
time_varying_covs = "age",
chains = 2, iter = 1000)
summary(fit)
}
}
\seealso{
\code{\link{ctbn_map}} for the fast L-BFGS-B / Laplace
approximation; \code{\link{ctbn_fit}} for the unified entry point.
}
