% Auto-generated --- do not edit by hand
\name{generate_true_interactions}
\alias{generate_true_interactions}
\title{Generate true interaction effects for an arbitrary order}
\usage{
generate_true_interactions(graph,
                           order = 2L,
                           prob = 0.4,
                           sd = 0.3,
                           strong_prob = 0.3,
                           seed = NULL)
}
\arguments{
\item{graph}{Named list of parent vectors (as in a network spec).}
\item{order}{Integer interaction order in [2, 5].}
\item{prob}{Probability of being non-zero.}
\item{sd}{Standard deviation of weak non-zero effects.}
\item{strong_prob}{Probability a non-zero effect is "strong".}
\item{seed}{Optional integer seed.}
}
\description{
Helper used internally by \code{\link{make_random_network}} but
exposed for users who want to extend an existing network with
higher-order interactions.
}
\value{
A named list of interaction coefficients keyed
\code{"<target>|<parent1>|<parent2>|..."}.
}
\examples{
net <- make_random_network(n_nodes = 6, max_parents = 3, seed = 1)
int3 <- generate_true_interactions(net$graph, order = 3,
                                    prob = 0.5, seed = 1)
str(int3, max.level = 1)

}
