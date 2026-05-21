% Auto-generated --- do not edit by hand
\name{make_default_network}
\alias{make_default_network}
\title{Default 10-node multimorbidity network}
\usage{
make_default_network(seed = 2026L)
}
\arguments{
\item{seed}{Integer, random seed for the stochastic 2-way and 3-way
interaction sampling step (to match the original simulation).}
}
\description{
Returns the canonical 10-condition network used in the original
CTBN multimorbidity studies (HTN, IHD, MN, OA, DM, HL, DVL, RS,
DE, CLRD), with hand-tuned baseline rates, parent sets, and effect
sizes drawn from UK Biobank empirical analyses.
}
\value{
A network specification list compatible with
\code{\link{simulate_ctbn_data}}.
}
\examples{
net <- make_default_network()
net$conds
}
