% Auto-generated --- do not edit by hand
\name{simulate_ctbn_data}
\alias{simulate_ctbn_data}
\title{Simulate a CTBN multimorbidity dataset}
\usage{
simulate_ctbn_data(network,
                   n_patients = 1000L,
                   t_horizon = 10,
                   dt_step = 5,
                   interaction_order = 3L,
                   use_smoking_dynamics = TRUE,
                   seed = NULL,
                   parallel = FALSE,
                   n_cores = 2L)
}
\arguments{
\item{network}{A network specification as returned by
\code{\link{make_default_network}} or \code{\link{make_random_network}}.}
\item{n_patients}{Integer, number of patients to simulate.}
\item{t_horizon}{Numeric, total follow-up duration (years).}
\item{dt_step}{Numeric, panel grid step (years).}
\item{interaction_order}{Integer in [1, 5]. Controls which interaction
terms enter the true rate function:
\itemize{
\item 1 -- main effects only
\item 2 -- main + 2-way
\item 3 -- main + 2-way + 3-way (default)
\item 4 -- up to 4-way
\item 5 -- up to 5-way
}}
\item{use_smoking_dynamics}{Logical, simulate the default smoking
transition process (passes as a time-varying covariate).}
\item{seed}{Optional integer random seed.}
\item{parallel}{Logical, parallelise patient simulation across cores.}
\item{n_cores}{Integer, number of cores when \code{parallel = TRUE}.}
}
\description{
Generates a panel-data dataset under a user-specified network. Each
patient's trajectory of conditions and time-varying covariates is
simulated forward in time using exact rate-based event sampling.

This generic generator works with any number of nodes (>= 2) and any
interaction order from 1 (main effects only) up to 5 (5-way
interactions), generalising the original 10-node multimorbidity
simulation. The user specifies the desired interaction order via the
\code{interaction_order} argument.
}
\value{
A long-format data.table with one row per at-risk interval
per patient. Columns: \code{eid}, \code{time_to_event}, \code{dt},
one column per condition, plus \code{age}, \code{sex_male},
\code{smk_current}, \code{smk_former}, and \code{event_cond}
(target of the next event, or NA for a censored interval).
}
\examples{
\dontrun{
# Use the default 10-node network with 2-way interactions
net <- make_default_network()
DT  <- simulate_ctbn_data(net, n_patients = 100,
interaction_order = 2, seed = 1)

# Random 6-node network with 3-way interactions
net2 <- make_random_network(n_nodes = 6, max_parents = 3, seed = 1)
DT2  <- simulate_ctbn_data(net2, n_patients = 200,
interaction_order = 3, seed = 1)
}
}
