% Auto-generated --- do not edit by hand
\name{print.ctbn_fit}
\alias{print.ctbn_fit}
\alias{print.ctbn_cph_fit}
\alias{print.ctbn_fctbn_fit}
\alias{print.ctbn_classic_fit}
\alias{print.ctbn_stan_fit}
\alias{print.ctbn_map_fit}
\title{Print a CTBN fit}
\usage{
\method{print}{ctbn_fit}(x, digits = 4, ...)
}
\arguments{
\item{x}{A \code{ctbn_fit} object.}
\item{digits}{Integer, decimal places.}
\item{...}{Currently unused.}
}
\description{
Print a CTBN fit
}
\value{
\code{x}, invisibly. Called for the side effect of printing.
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 1)

print(fit)
}

}
