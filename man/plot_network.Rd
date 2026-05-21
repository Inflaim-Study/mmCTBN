% Auto-generated --- do not edit by hand
\name{plot_network}
\alias{plot_network}
\title{Plot a fitted CTBN as a directed network}
\usage{
plot_network(fit,
             pip_threshold = NULL,
             layout = "fr",
             node_size = 12,
             show_intensity_labels = TRUE,
             intensity_scale = 1000,
             seed = 123)
}
\arguments{
\item{fit}{A \code{ctbn_fit} object.}
\item{pip_threshold}{Selection threshold (default
\code{fit$call_args$pip_threshold} or 0.5).}
\item{layout}{Character; igraph/ggraph layout name (default \code{"fr"},
the Fruchterman-Reingold force-directed layout).}
\item{node_size}{Numeric, node size (default 12).}
\item{show_intensity_labels}{Logical; if TRUE (default), edges are
labelled with their estimated transition intensity (per
\code{intensity_scale} units of person-time).}
\item{intensity_scale}{Numeric, the amount of person-time the displayed
edge intensities are scaled to (default 1000).}
\item{seed}{Integer; random seed passed to \code{set.seed()} so the
force-directed layout is reproducible (default 123).}
}
\description{
Draws a DAG of significant influencer -> target edges. Edge colour
encodes log rate-ratio (red = excess risk, blue = protective); edge
width encodes |log RR|.
}
\value{
A \code{ggplot} object (if \pkg{ggraph} available) or
invisibly the active-edge data.frame.
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 2)

# Directed network graph (requires ggraph + igraph).
plot_network(fit, layout = "circle")
}

}
