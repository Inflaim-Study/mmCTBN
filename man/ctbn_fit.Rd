% Auto-generated --- do not edit by hand
\name{ctbn_fit}
\alias{ctbn_fit}
\title{Fit a Continuous-Time Bayesian Network (unified entry point)}
\usage{
ctbn_fit(DT_wide,
         method = c("map", "stan", "classic", "fctbn", "cph"),
         prior = c("spike_slab", "structured", "lasso", "horseshoe"),
         max_order = 1L,
         fixed_covs = character(0),
         time_varying_covs = character(0),
         target_conditions = NULL,
         variable_select = FALSE,
         pip_threshold = 0.5,
         theta = 1.0,
         parallel = FALSE,
         verbose = TRUE,
         ...)
}
\arguments{
\item{DT_wide}{A data.table in wide format. Must contain columns
\code{eid}, \code{time_to_event}, \code{dt}, plus one column per
condition and per covariate. See \code{\link{simulate_ctbn_data}}
for an example of the expected layout.}
\item{method}{Character; one of \code{"map"} (default),
\code{"stan"}, \code{"classic"}, \code{"fctbn"}, \code{"cph"}.}
\item{prior}{Prior family used by \code{method="map"} and
\code{method="stan"}: \code{"spike_slab"} (default),
\code{"structured"}, \code{"lasso"}, or \code{"horseshoe"}.
Ignored by the other methods.}
\item{max_order}{Integer in [1, 5]. Maximum interaction order for
the MAP/Stan engines. The classical, FCTBN and CPH methods always
use \code{max_order = 1} (main effects only on the parent side).}
\item{fixed_covs}{Character vector of fixed-covariate column names.}
\item{time_varying_covs}{Character vector of time-varying covariate
column names.}
\item{target_conditions}{Optional subset of conditions to model as
outcomes.}
\item{variable_select}{Used by MAP/Stan only;
silently ignored otherwise.}
\item{pip_threshold}{Used by MAP/Stan only;
silently ignored otherwise.}
\item{theta}{Used by MAP/Stan only;
silently ignored otherwise.}
\item{parallel}{Common across methods.}
\item{verbose}{Common across methods.}
\item{...}{Additional method-specific arguments. Examples:
\itemize{
\item \code{method = "stan"}: \code{chains}, \code{iter},
\code{warmup}, \code{seed}.
\item \code{method = "map"}: \code{lbfgs_maxit},
\code{compute_se}, \code{a0}, \code{b0}, \dots
\item \code{method = "classic"} or \code{"cph"}: \code{parents},
\code{max_parents}, \code{add_laplace}; \code{"cph"} also
accepts \code{alpha} (p-value gate) and \code{ties}.
\item \code{method = "fctbn"}: \code{lambda}, \code{max_iter},
\code{tol}, \code{gmm_stop}, \code{gmm_warmup}, \code{adaptive},
\code{pilot_iter}.
}}
}
\description{
Single-call wrapper that dispatches to one of five estimation
backends:
\describe{
\item{\code{"map"}}{Regularised Poisson MAP + Laplace SE
(\code{\link{ctbn_map}}). Fast, supports up to 5-way interactions
and four prior families.}
\item{\code{"stan"}}{Full Bayesian MCMC (\code{\link{ctbn_stan}}).
Slower but provides full posteriors.}
\item{\code{"classic"}}{Classical CTBN of Nodelman et al.
(\code{\link{ctbn_classic}}). Fully parametric, no covariates;
closed-form MLE per parent configuration.}
\item{\code{"fctbn"}}{Functional CTBN of Faruqui et al.
(\code{\link{ctbn_fctbn}}). Poisson regression with adaptive
group lasso and multiplicative parent effects.}
\item{\code{"cph"}}{CTBN-CPH of Guillamet et al.
(\code{\link{ctbn_cph}}). Classical baseline rates modulated by
per-edge Cox proportional-hazards regressions, yielding
individualised intensities.}
}

All backends produce an object that inherits from class
\code{"ctbn_fit"}, so downstream tools (\code{\link{get_lp}},
\code{\link{ctbn_cv}}, \code{\link{compute_F_m}},
\code{\link{compute_interaction_effects}}, the plotting suite, and
the metric helpers) work uniformly across methods.
}
\value{
An object inheriting from \code{"ctbn_fit"}, with subclass
determined by \code{method}.
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500,
interaction_order = 2, seed = 1)

fit_map     <- ctbn_fit(DT, method = "map",
prior = "spike_slab", max_order = 2,
fixed_covs = c("age", "sex_male"))
fit_classic <- ctbn_fit(DT, method = "classic", max_parents = 3)
fit_fctbn   <- ctbn_fit(DT, method = "fctbn",
fixed_covs = c("age", "sex_male"),
max_iter = 5000)
fit_cph     <- ctbn_fit(DT, method = "cph",
fixed_covs = c("age", "sex_male"),
max_parents = 3, alpha = 0.05)
}
}
\seealso{
\code{\link{ctbn_map}}, \code{\link{ctbn_stan}},
\code{\link{ctbn_classic}}, \code{\link{ctbn_fctbn}},
\code{\link{ctbn_cph}}, \code{\link{simulate_ctbn_data}},
\code{\link{ctbn_cv}}.
}
