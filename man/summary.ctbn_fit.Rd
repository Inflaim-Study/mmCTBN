% Auto-generated --- do not edit by hand
\name{summary.ctbn_fit}
\alias{summary.ctbn_fit}
\alias{summary.ctbn_cph_fit}
\alias{summary.ctbn_fctbn_fit}
\alias{summary.ctbn_classic_fit}
\alias{summary.ctbn_stan_fit}
\alias{summary.ctbn_map_fit}
\title{Summary of active influencer pairs}
\usage{
\method{summary}{ctbn_fit}(object,
                           pip_threshold = NULL,
                           show_covariates = FALSE,
                           ...)
}
\arguments{
\item{object}{A \code{ctbn_fit} object.}
\item{pip_threshold}{Selection threshold; default uses
\code{object$call_args$pip_threshold}.}
\item{show_covariates}{Logical, also list active covariates.}
\item{...}{Currently unused.}
}
\description{
Summary of active influencer pairs
}
\value{
A list summarising the active influencer pairs, returned
invisibly; printed for inspection.
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 1)

summary(fit)
}

}
