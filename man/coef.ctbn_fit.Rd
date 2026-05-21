% Auto-generated --- do not edit by hand
\name{coef.ctbn_fit}
\alias{coef.ctbn_fit}
\title{Extract coefficients}
\usage{
\method{coef}{ctbn_fit}(object,
                        type = c("beta", "intensity", "pip", "kappa", "se"),
                        ...)
}
\arguments{
\item{object}{A \code{ctbn_fit} object.}
\item{type}{Which coefficient table to return:
\code{"beta"} (default; condition-condition log-RR),
\code{"intensity"}, \code{"pip"}, \code{"kappa"}, or \code{"se"}.}
\item{...}{Currently unused.}
}
\description{
Extract coefficients
}
\value{
A numeric matrix of the requested coefficient table (may be
\code{NULL} if the fit does not provide it).
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 1)

coef(fit)                 # condition-condition log-RR matrix
coef(fit, "intensity")    # reference-profile intensities
}

}
