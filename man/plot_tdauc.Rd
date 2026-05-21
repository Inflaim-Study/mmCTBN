% Auto-generated --- do not edit by hand
\name{plot_tdauc}
\alias{plot_tdauc}
\title{Plot time-dependent AUC from cross-validation}
\usage{
plot_tdauc(cv_results, targets = NULL)
}
\arguments{
\item{cv_results}{A data.table from \code{\link{ctbn_cv}}.}
\item{targets}{Optional character subset to display.}
}
\description{
Faceted line plot of mean tdAUC vs eval_time, with a 95% across-fold
ribbon, for each (model, target).
}
\value{
A \code{ggplot} object.
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 400, seed = 1)
DTw <- prepare_wide(DT)

fit_fns <- list(map = function(DT, ...) ctbn_map(DT, max_order = 1, ...))
cv      <- ctbn_cv(DTw, fit_fns = fit_fns, k_folds = 3,
                    eval_times = c(1, 3, 5), seed = 1)

plot_tdauc(cv)
}

}
