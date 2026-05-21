% Auto-generated --- do not edit by hand
\name{make_random_network}
\alias{make_random_network}
\title{Generate a random network specification for simulation}
\usage{
make_random_network(n_nodes = 10L,
                    min_parents = 1L,
                    max_parents = 4L,
                    node_names = NULL,
                    baseline_range = c(0.005, 0.20),
                    main_effect_sd = 0.35,
                    int2_prob = 0.40,
                    int3_prob = 0.30,
                    int4_prob = 0.20,
                    int5_prob = 0.10,
                    int_effect_sd = 0.30,
                    strong_int_prob = 0.30,
                    covariate_names = c("age", "sex_male", "smk_current", "smk_former"),
                    seed = NULL)
}
\arguments{
\item{n_nodes}{Integer, number of conditions in the network (>= 2).}
\item{min_parents}{Integer, minimum number of parents per condition.}
\item{max_parents}{Integer, maximum number of parents per condition.
Must be <= n_nodes - 1.}
\item{node_names}{Character vector of length n_nodes giving custom
node names. Defaults to \code{C1, C2, ...}.}
\item{baseline_range}{Numeric length-2 vector; baseline rates are
drawn uniformly from this range (then log-transformed).}
\item{main_effect_sd}{Numeric, standard deviation of main-effect betas.}
\item{int2_prob}{Numeric in [0,1];
marginal probability that any 2-way / 3-way / 4-way / 5-way
interaction term is non-zero.}
\item{int3_prob}{Numeric in [0,1];
marginal probability that any 2-way / 3-way / 4-way / 5-way
interaction term is non-zero.}
\item{int4_prob}{Numeric in [0,1];
marginal probability that any 2-way / 3-way / 4-way / 5-way
interaction term is non-zero.}
\item{int5_prob}{Numeric in [0,1];
marginal probability that any 2-way / 3-way / 4-way / 5-way
interaction term is non-zero.}
\item{int_effect_sd}{Numeric, SD of non-zero interaction betas.}
\item{strong_int_prob}{Probability that a non-zero interaction is a
"strong" effect (mean ~1.0) rather than a "weak" one (mean 0).}
\item{covariate_names}{Character vector of covariate names (passed
through to gamma generation).}
\item{seed}{Integer random seed.}
}
\description{
Creates a network of `n_nodes` conditions, where each condition has
between `min_parents` and `max_parents` parent influencers chosen at
random. All effect sizes (baseline rates, main effects, 2-way to
5-way interactions, covariate gammas) are sampled from sensible
defaults so the resulting simulated data mimic real EHR multimorbidity
patterns.
}
\value{
A list with components:
\describe{
\item{conds}{Character vector of condition names.}
\item{beta0}{Named numeric vector of log baseline rates.}
\item{graph}{Named list: parents per condition.}
\item{beta_main}{Named numeric vector of main-effect coefficients.}
\item{beta_int2}{Named list of 2-way interaction coefficients.}
\item{beta_int3}{Named list of 3-way coefficients (if any).}
\item{beta_int4}{Named list of 4-way coefficients (if any).}
\item{beta_int5}{Named list of 5-way coefficients (if any).}
\item{gamma}{Named list of per-node covariate coefficients.}
\item{init_prev}{Named numeric vector of initial prevalences.}
}
}
\examples{
net <- make_random_network(n_nodes = 8, max_parents = 3, seed = 1)
str(net, max.level = 1)
}
