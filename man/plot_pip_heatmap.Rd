% Auto-generated --- do not edit by hand
\name{plot_pip_heatmap}
\alias{plot_pip_heatmap}
\title{Heatmap of posterior inclusion / pseudo-PIP / shrinkage}
\usage{
plot_pip_heatmap(fit, threshold = 0.5, show_values = TRUE)
}
\arguments{
\item{fit}{A \code{ctbn_fit} object.}
\item{threshold}{Numeric in [0, 1]; inclusion-score cut-off (default 0.5).
Cells whose score is at or above this value are highlighted with bold
white labels.}
\item{show_values}{Logical; if TRUE (default), numeric scores at or above
\code{threshold} are overlaid on the heatmap.}
}
\description{
For \code{spike_slab} and \code{structured} priors, visualises the
inclusion probability matrix; for \code{lasso}/\code{horseshoe},
visualises \code{1 - kappa} (a shrinkage-based inclusion proxy).
}
\value{
A \code{ggplot} object.
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 1,
                 variable_select = TRUE)

plot_pip_heatmap(fit, threshold = 0.5)
}

}
