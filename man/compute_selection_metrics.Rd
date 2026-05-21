% Auto-generated --- do not edit by hand
\name{compute_selection_metrics}
\alias{compute_selection_metrics}
\title{Selection metrics (TPR / FPR / AUC for variable selection)}
\usage{
compute_selection_metrics(fit,
                          network,
                          pip_thresh = 0.5,
                          max_order_truth = NULL)
}
\arguments{
\item{fit}{A fitted \code{ctbn_fit} with selection scores.}
\item{network}{The true network specification.}
\item{pip_thresh}{Selection threshold (default 0.5).}
\item{max_order_truth}{Highest-order interactions to evaluate
(>= fit's max_order is fine).}
}
\description{
Selection metrics (TPR / FPR / AUC for variable selection)
}
\value{
A data.table with one row per interaction order
(\code{main_effect}, \code{interaction_2way}, \code{interaction_3way},
...) and columns \code{tpr}, \code{fpr}, \code{selection_auc}.
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500,
                           interaction_order = 2, seed = 1)
DTw <- prepare_wide(DT)

fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 2,
                 variable_select = TRUE)
compute_selection_metrics(fit, network = net, pip_thresh = 0.5)
}

}
