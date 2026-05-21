% Auto-generated --- do not edit by hand
\name{ctbn_cph}
\alias{ctbn_cph}
\title{Fit a CTBN-CPH (Guillamet 2025) — individualised CTBN via Cox-PH}
\usage{
ctbn_cph(DT_wide,
         parents = NULL,
         max_parents = 4L,
         fixed_covs = character(0),
         time_varying_covs = character(0),
         target_conditions = NULL,
         alpha = 0.05,
         add_laplace = 0.5,
         ties = "efron",
         parallel = FALSE,
         verbose = TRUE)
}
\arguments{
\item{DT_wide}{A wide-format data.table; same layout as
\code{\link{ctbn_map}}.}
\item{parents}{Named list of parent vectors per condition (the
network topology). If \code{NULL}, every other condition is
treated as a parent.}
\item{max_parents}{Integer cap on per-target parent set size.}
\item{fixed_covs}{Covariate columns to enter the
Cox stage.}
\item{time_varying_covs}{Covariate columns to enter the
Cox stage.}
\item{target_conditions}{Optional subset of conditions to model.}
\item{alpha}{Numeric significance threshold on Cox p-values; any
covariate with \eqn{p \ge \alpha} is gated out at inference.}
\item{add_laplace}{Pseudo-count for the classical baseline stage.}
\item{ties}{Character; ties option passed to
\code{survival::coxph}. Default \code{"efron"}.}
\item{parallel}{As elsewhere.}
\item{verbose}{As elsewhere.}
}
\description{
Hybrid two-stage model. First, baseline rates
\eqn{\lambda^0_{m\,|\,u}} are estimated with a classical CTBN (see
\code{\link{ctbn_classic}}). Second, for every directed edge
\eqn{j \to m} in the topology a Cox proportional-hazards regression
is run on the time-to-transition of \eqn{m}, conditioning on
patient covariates. At inference time the intensity for a patient
with covariates \eqn{C_i} is
\deqn{\lambda(C_i;\, m\,|\,u) \;=\; \lambda^0_{m\,|\,u}
\cdot \exp\!\Big( \sum_{l} \beta_l\, c_l\, \sigma_l \Big),}
where \eqn{\sigma_l = 1} iff the Cox p-value for covariate \eqn{l}
on edge \eqn{j \to m} is below \code{alpha}.
}
\value{
An object of class \code{c("ctbn_cph_fit", "ctbn_fit")}
with the standard slots plus
\describe{
\item{classic}{The underlying classical CTBN fit (baseline
lambdas).}
\item{cox_fits}{Nested list: \code{cox_fits[[target]][[parent]]}
is a Cox model fitted to time-to-target conditional on parent
presence.}
\item{beta_cov}{Nested list of named numeric vectors of Cox
coefficients per edge.}
\item{p_cov}{Nested list of named numeric vectors of Cox
p-values per edge.}
}
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
DTw <- prepare_wide(DT)

# CTBN-CPH: classical baseline modulated by per-edge Cox regressions.
# Requires the 'survival' package.
fit <- ctbn_cph(DTw, fixed_covs = c("age", "sex_male"),
                 max_parents = 3, alpha = 0.05)
summary(fit)
}

}
