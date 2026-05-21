% Auto-generated --- do not edit by hand
\name{build_design_matrix}
\alias{build_design_matrix}
\title{Build the per-node design matrix Phi_m}
\usage{
build_design_matrix(X_all, Z, cond_names, target, max_order)
}
\arguments{
\item{X_all}{Numeric matrix of all condition indicators
(rows = at-risk intervals, columns = condition names).}
\item{Z}{Numeric matrix of covariates
(rows = at-risk intervals, columns = covariate names).}
\item{cond_names}{Character vector of all condition names.}
\item{target}{Character, name of the target condition.}
\item{max_order}{Integer, maximum interaction order
(1 = main effects only; up to 5 allowed).}
}
\description{
Assembles the full design matrix
\eqn{\Phi_m = [X_{-m} \mid \text{interactions} \mid Z]}
for a single target node \eqn{m}, ready for use in MAP/MCMC fitting.
}
\value{
A list with components:
\describe{
\item{Phi}{Numeric matrix; full design matrix for target.}
\item{x_cols}{Character vector of column names for the X-block.}
\item{x_orders}{Integer vector of interaction orders for X columns
(0 = main effect, 1 = 2-way, 2 = 3-way, ...).}
\item{n_beta}{Number of X columns (condition + interaction terms).}
\item{n_gamma}{Number of Z columns (covariates).}
}
}
\examples{
# Two conditions plus one covariate, target = "B".
X <- matrix(c(1, 0, 1,
              0, 1, 1), ncol = 2,
            dimnames = list(NULL, c("A", "B")))
Z <- matrix(c(0.2, -0.5, 1.1), ncol = 1, dimnames = list(NULL, "age"))
dm <- build_design_matrix(X, Z, cond_names = c("A", "B"),
                           target = "B", max_order = 1)
dm$x_cols
dm$Phi

}
