% Auto-generated --- do not edit by hand
\name{prepare_wide}
\alias{prepare_wide}
\title{Prepare a long-format CTBN dataset into wide format}
\usage{
prepare_wide(DT, cond_names = NULL)
}
\arguments{
\item{DT}{A data.table in long (interval) format produced by
\code{\link{simulate_ctbn_data}} or compiled from EHR data.}
\item{cond_names}{Character vector of condition column names. If
\code{NULL}, all integer/0-1 columns excluding the reserved fields
are taken to be conditions.}
}
\description{
Strips the \code{event_cond} bookkeeping column (used during data
generation) and ensures all condition columns are integers. The
result is suitable for input to \code{\link{ctbn_fit}}.
}
\value{
A data.table in the canonical wide format expected by
\code{ctbn_fit}.
}
\examples{
net <- make_random_network(n_nodes = 5, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 50, seed = 1)
DTw <- prepare_wide(DT)
head(DTw)

}
