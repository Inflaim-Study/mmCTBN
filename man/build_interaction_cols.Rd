% Auto-generated --- do not edit by hand
\name{build_interaction_cols}
\alias{build_interaction_cols}
\title{Add interaction columns to a data.frame}
\usage{
build_interaction_cols(data, vars, k)
}
\arguments{
\item{data}{A data.frame or data.table containing the base variables.}
\item{vars}{Character vector of base variable names to interact.}
\item{k}{Integer, maximum interaction order (1 = main effects only,
2 = include 2-way, ..., up to length(vars)).}
}
\description{
Constructs all pairwise, three-way, ..., k-way interaction columns
from a set of base variables and appends them to the data.frame
(or data.table). Each new column is named by joining the
parent variable names with `"_x_"`.
}
\value{
A list with elements:
\describe{
\item{data}{The augmented data.frame including the new columns.}
\item{cols}{Character vector of new interaction column names.}
\item{orders}{Integer vector of interaction orders (1 = 2-way,
2 = 3-way, ...).}
}
}
\examples{
df <- data.frame(A = c(1, 0, 1), B = c(1, 1, 0), C = c(0, 1, 1))
res <- build_interaction_cols(df, c("A", "B", "C"), k = 3)
res$cols
}
