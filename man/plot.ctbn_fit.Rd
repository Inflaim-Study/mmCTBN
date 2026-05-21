% Auto-generated --- do not edit by hand
\name{plot.ctbn_fit}
\alias{plot.ctbn_fit}
\title{Plot a fitted CTBN}
\usage{
\method{plot}{ctbn_fit}(x, type = c("network","heatmap","pip"), ...)
}
\arguments{
\item{x}{A \code{ctbn_fit} object.}
\item{type}{One of \code{"network"} (DAG), \code{"heatmap"} (RR
heatmap), \code{"pip"} (PIP heatmap).}
\item{...}{Passed to the underlying plot function.}
}
\description{
Default \code{plot} dispatches to \code{\link{plot_network}}.
}
\value{
A \code{ggplot} object from the underlying plot function.
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 1)

plot(fit, type = "network")
plot(fit, type = "heatmap")
}

}
