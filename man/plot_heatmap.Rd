% Auto-generated --- do not edit by hand
\name{plot_heatmap}
\alias{plot_heatmap}
\title{Heatmap of pairwise rate ratios}
\usage{
plot_heatmap(fit,
             scale = c("RR", "logRR"),
             mask_inactive = TRUE,
             pip_threshold = NULL,
             show_values = TRUE)
}
\arguments{
\item{fit}{A \code{ctbn_fit} object.}
\item{scale}{One of \code{"RR"} (default) or \code{"logRR"}.}
\item{mask_inactive}{Logical; if TRUE (default), cells that fail the
selection rule are greyed out.}
\item{pip_threshold}{Selection threshold (default
\code{fit$call_args$pip_threshold} or 0.5).}
\item{show_values}{Logical; if TRUE (default), the numeric value of each
active cell is overlaid on the heatmap.}
}
\description{
Visualises the \code{beta_matrix} (log-RR) of a fitted CTBN as a
diverging heatmap. Diagonals are blanked.
}
\value{
A \code{ggplot} object.
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 1)

plot_heatmap(fit, scale = "RR")
}

}
