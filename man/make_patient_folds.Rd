% Auto-generated --- do not edit by hand
\name{make_patient_folds}
\alias{make_patient_folds}
\title{Construct patient-level fold assignments}
\usage{
make_patient_folds(DT, k = 5, seed = 42)
}
\arguments{
\item{DT}{A data.table with an \code{eid} column.}
\item{k}{Integer, number of folds (default 5).}
\item{seed}{Integer, RNG seed (default 42).}
}
\description{
Splits patients (\code{eid}) — not rows — into \code{k} folds. This
avoids leakage of within-patient correlation between train and test.
}
\value{
A data.table with columns \code{eid} and \code{fold}.
}
\examples{
net <- make_random_network(n_nodes = 5, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 60, seed = 1)
folds <- make_patient_folds(DT, k = 3, seed = 1)
head(folds)
table(folds$fold)

}
