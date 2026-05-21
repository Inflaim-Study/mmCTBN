% Auto-generated --- do not edit by hand
\name{compute_F_m}
\alias{compute_F_m}
\title{Cumulative incidence F_m(tau) for a target node and patient profile}
\usage{
compute_F_m(fit, target, patient_dt, tau, zbar = NULL, unit = 1)
}
\arguments{
\item{fit}{A fitted \code{ctbn_fit} object.}
\item{target}{Character, the target condition.}
\item{patient_dt}{A data.table of at-risk intervals for one patient.}
\item{tau}{Numeric, time horizon at which to evaluate F_m.}
\item{zbar}{Optional named numeric vector to override covariate values.}
\item{unit}{Optional rate scaling factor.}
}
\description{
Computes \eqn{F_m(\tau) = 1 - \exp(-\int_0^\tau \lambda_m(s)\, ds)}
for a given patient trajectory, evaluated under a fitted CTBN. The
rate is plugged in at the posterior mean (Stan) or MAP estimate
(MAP fit). Works with both \code{ctbn_map_fit} and \code{ctbn_stan_fit}.
}
\value{
Numeric, cumulative incidence in [0, 1].
}
\examples{
\dontrun{
net <- make_random_network(n_nodes = 6, seed = 1)
DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
DTw <- prepare_wide(DT)
fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 2)

# Cumulative incidence of condition "C1" by 5 years for one patient.
target  <- fit$call_args$all_conditions[1]
pat_dt  <- DTw[eid == DTw$eid[1]]
compute_F_m(fit, target = target, patient_dt = pat_dt, tau = 5)
}

}
