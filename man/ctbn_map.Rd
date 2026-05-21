% Auto-generated --- do not edit by hand
\name{ctbn_map}
\alias{ctbn_map}
\title{Fit a CTBN by Maximum-A-Posteriori (MAP) estimation}
\usage{
ctbn_map(DT_wide,
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
         sigma2_init = NULL,
         lambda2_init = NULL,
         sigma2_g_init = NULL,
         lbfgs_maxit = 500L,
         compute_se = TRUE,
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
\item{a0}{Hyperparameters of the inverse-gamma slab prior on
\eqn{\sigma^2} (structured / lasso / spike_slab).}
\item{b0}{Hyperparameters of the inverse-gamma slab prior on
\eqn{\sigma^2} (structured / lasso / spike_slab).}
\item{pi0}{Base inclusion probability (spike_slab only).}
\item{spike_var}{Spike variance (spike_slab only).}
\item{tau0}{Half-Cauchy global scale (horseshoe only).}
\item{slab_df}{Slab degrees of freedom (horseshoe only).}
\item{slab_scale}{Slab scale (horseshoe only).}
\item{sigma2_init}{Optional manual overrides
for the empirical-Bayes hyperparameter initialisation.}
\item{lambda2_init}{Optional manual overrides
for the empirical-Bayes hyperparameter initialisation.}
\item{sigma2_g_init}{Optional manual overrides
for the empirical-Bayes hyperparameter initialisation.}
\item{lbfgs_maxit}{Integer, maximum L-BFGS-B iterations.}
\item{compute_se}{Logical; compute Laplace SEs (default TRUE).}
\item{parallel}{Logical; if TRUE distribute target nodes across workers
via \pkg{future.apply}. Requires \code{plan(multisession)} or similar
to be set by the caller.}
\item{verbose}{Logical, print progress messages.}
}
\description{
Fast L-BFGS-B optimiser for the joint posterior mode plus a Laplace
approximation for posterior standard errors. Provides a 50-200x speed-up
over the equivalent Stan MCMC fit (\code{\link{ctbn_stan}}) at the cost
of replacing posterior samples with normal approximations centred at
the MAP.
}
\value{
An object of class \code{c("ctbn_map_fit", "ctbn_fit")} with
slots \code{beta_matrix}, \code{se_matrix}, \code{pip_matrix},
\code{kappa_matrix}, \code{intensity_matrix}, \code{pip_cov_list},
\code{convergence}, \code{map_fits}, \code{call_args}.
}
\examples{
\dontrun{
# Simulate from a 6-node network with 2-way interactions
net <- make_random_network(n_nodes = 6, max_parents = 2, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500,
interaction_order = 2, seed = 1)
DTw <- prepare_wide(DT)

# Fast MAP fit
fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 2,
fixed_covs = c("intercept","sex_male"),
time_varying_covs = "age")
summary(fit)
}
}
\seealso{
\code{\link{ctbn_fit}} for the unified entry point;
\code{\link{ctbn_stan}} for full Bayesian MCMC.
}
