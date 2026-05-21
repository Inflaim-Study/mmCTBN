% Auto-generated --- do not edit by hand
\name{get_lp}
\alias{get_lp}
\alias{get_lp.ctbn_fit}
\alias{get_lp.ctbn_stan_fit}
\alias{get_lp.ctbn_map_fit}
\alias{get_lp.ctbn_fctbn_fit}
\alias{get_lp.ctbn_cph_fit}
\alias{get_lp.ctbn_classic_fit}
\title{Linear predictor and rate for a CTBN fit}
\usage{
get_lp(fit, newdata, target, ...)
}
\arguments{
\item{fit}{A fitted \code{ctbn_fit} object.}
\item{newdata}{A data.table of at-risk intervals.}
\item{target}{Character, the target condition.}
\item{...}{Currently unused.}
}
\description{
Generic that returns the linear predictor \eqn{\eta = X\beta + Z\gamma}
and rate \eqn{\lambda = \exp(\eta)} for one target node, evaluated on
new data. S3 methods dispatch on whether \code{fit} is a MAP fit
(\code{ctbn_map_fit}) or a Stan fit (\code{ctbn_stan_fit}).
}
\value{
A list with elements \code{lp} (linear predictor) and
\code{lambda} (rate).
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
DTw <- prepare_wide(DT)
fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 1)

target <- fit$call_args$all_conditions[1]
pred   <- get_lp(fit, DTw, target)
str(pred)
}

}
