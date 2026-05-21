# =============================================================================
# priors.R  --  Prior log-density, gradient and Hessian functions
# =============================================================================
#
# Four prior families are supported:
#   1. structured  -- Order-penalised normal:   beta_j ~ N(0, sigma^2_p * exp(-theta*p))
#   2. lasso       -- Smooth Laplace:           beta_j ~ DE(0, s_p)
#   3. spike_slab  -- Marginalised mixture:     beta_j | delta_j ~ N0/N_slab
#   4. horseshoe   -- Regularised MAP approx.   beta_j ~ N(0, sigma^2_hs(p))
#
# Each family exposes three internal helpers:
#   .prior_logdens_<name>()  log p(beta | hyperparams)
#   .grad_prior_<name>()     gradient w.r.t. beta
#   .hess_diag_prior_<name>() Hessian diagonal w.r.t. beta
#
# The order penalty parameter theta scales how aggressively higher-order
# interactions are shrunk relative to main effects.
# =============================================================================


# ── 1. Structured (Order-penalised Normal) ───────────────────────────────────

.prior_logdens_structured <- function(beta, sigma2, orders, theta) {
  v <- sigma2[orders + 1L] * exp(-theta * orders)
  -0.5 * sum(beta^2 / v + log(2 * pi * v))
}

.grad_prior_structured <- function(beta, sigma2, orders, theta) {
  v <- sigma2[orders + 1L] * exp(-theta * orders)
  -beta / v
}

.hess_diag_prior_structured <- function(sigma2, orders, theta) {
  v <- sigma2[orders + 1L] * exp(-theta * orders)
  -1.0 / v
}


# ── 2. LASSO (Smooth Laplace) ────────────────────────────────────────────────

.prior_logdens_lasso <- function(beta, lambda2, orders, theta, eps = 1e-6) {
  s <- sqrt(lambda2[orders + 1L]) * exp(-theta * orders / 2.0)
  -sum(sqrt(beta^2 + eps^2) / s + log(2 * s))
}

.grad_prior_lasso <- function(beta, lambda2, orders, theta, eps = 1e-6) {
  s <- sqrt(lambda2[orders + 1L]) * exp(-theta * orders / 2.0)
  -beta / (s * sqrt(beta^2 + eps^2))
}

.hess_diag_prior_lasso <- function(beta, lambda2, orders, theta, eps = 1e-6) {
  s <- sqrt(lambda2[orders + 1L]) * exp(-theta * orders / 2.0)
  -(eps^2) / (s * (beta^2 + eps^2)^1.5)
}


# ── 3. Spike-and-Slab (Marginalised Mixture) ─────────────────────────────────

.prior_logdens_spikeslab <- function(beta, sigma2, orders, theta,
                                     pi0, spike_var) {
  pi_j   <- pmin(pmax(pi0 * exp(-theta * orders), 1e-10), 1 - 1e-10)
  slab_v <- sigma2[orders + 1L]
  log_s  <- log(pi_j)      + stats::dnorm(beta, 0, sqrt(slab_v),    log = TRUE)
  log_k  <- log(1 - pi_j)  + stats::dnorm(beta, 0, sqrt(spike_var), log = TRUE)
  sum(apply(cbind(log_s, log_k), 1, function(r) {
    mx <- max(r); mx + log(sum(exp(r - mx)))
  }))
}

.grad_prior_spikeslab <- function(beta, sigma2, orders, theta, pi0, spike_var) {
  pi_j   <- pmin(pmax(pi0 * exp(-theta * orders), 1e-10), 1 - 1e-10)
  slab_v <- sigma2[orders + 1L]
  phi_s  <- stats::dnorm(beta, 0, sqrt(slab_v))
  phi_k  <- stats::dnorm(beta, 0, sqrt(spike_var))
  numer  <- -pi_j * phi_s * beta / slab_v -
             (1 - pi_j) * phi_k * beta / spike_var
  denom  <- pi_j * phi_s + (1 - pi_j) * phi_k + 1e-300
  numer / denom
}

.hess_diag_prior_spikeslab <- function(beta, sigma2, orders, theta,
                                       pi0, spike_var, eps = 1e-5) {
  gp <- .grad_prior_spikeslab(beta + eps, sigma2, orders, theta, pi0, spike_var)
  gm <- .grad_prior_spikeslab(beta - eps, sigma2, orders, theta, pi0, spike_var)
  (gp - gm) / (2 * eps)
}


# ── 4. Horseshoe (Regularised; MAP approximation) ────────────────────────────

.hs_eff_var <- function(hs_tau, hs_c2, orders, theta) {
  tau_j <- hs_tau * exp(-theta * orders)
  hs_c2 * tau_j^2 / (hs_c2 + tau_j^2)
}

.prior_logdens_horseshoe <- function(beta, hs_tau, hs_c2, orders, theta) {
  v <- .hs_eff_var(hs_tau, hs_c2, orders, theta)
  -0.5 * sum(beta^2 / v + log(2 * pi * v))
}

.grad_prior_horseshoe <- function(beta, hs_tau, hs_c2, orders, theta) {
  v <- .hs_eff_var(hs_tau, hs_c2, orders, theta)
  -beta / v
}

.hess_diag_prior_horseshoe <- function(hs_tau, hs_c2, orders, theta) {
  v <- .hs_eff_var(hs_tau, hs_c2, orders, theta)
  -1.0 / v
}

#' Compute E[c^2] under InvGamma(slab_df/2, slab_df * slab_scale^2 / 2)
#'
#' Used to compute the effective slab variance for the horseshoe MAP
#' approximation. Mean exists for slab_df > 2; default slab_df = 4
#' gives E[c^2] = 2 * slab_scale^2.
#'
#' @keywords internal
#' @noRd
.hs_c2_mean <- function(slab_df, slab_scale) {
  if (slab_df > 2)
    (slab_df / 2 * slab_scale^2) / (slab_df / 2 - 1)
  else
    slab_scale^2 * 4
}


# =============================================================================
# Poisson log-likelihood, gradient, Hessian (for at-risk panel data)
# =============================================================================

.poisson_loglik <- function(par, Phi, y, offset) {
  eta <- pmin(drop(Phi %*% par) + offset, 30)
  sum(y * eta - exp(eta))
}

.poisson_grad <- function(par, Phi, y, offset) {
  eta <- pmin(drop(Phi %*% par) + offset, 30)
  drop(crossprod(Phi, y - exp(eta)))
}

# Returns Phi^T diag(mu) Phi (Fisher information; PSD)
.poisson_neg_hess <- function(par, Phi, offset) {
  eta <- pmin(drop(Phi %*% par) + offset, 30)
  crossprod(Phi * sqrt(exp(eta)))
}


# =============================================================================
# Empirical Bayes initialisation of hyperparameters
# =============================================================================

.init_hyperparams_eb <- function(DT_wide, all_conds, all_covs,
                                  max_order, a0 = 2, b0 = 1) {
  betas_crude <- numeric(0)
  for (target in all_conds) {
    infl <- setdiff(all_conds, target)
    atr  <- DT_wide[get(target) == 0 & dt > 0]
    if (nrow(atr) < 10) next
    df  <- as.data.frame(atr)
    ec  <- paste0(target, "_event")
    if (!ec %in% names(df)) next
    X_  <- as.matrix(df[, infl, drop = FALSE])
    y_  <- as.integer(df[[ec]])
    off <- log(pmax(df$dt, 1e-10))
    tryCatch({
      fit_ <- stats::glm.fit(cbind(1, X_), y_,
                              family  = stats::poisson(link = "log"),
                              offset  = off,
                              control = list(maxit = 30))
      betas_crude <- c(betas_crude, stats::coef(fit_)[-1])
    }, error = function(e) NULL)
  }
  emp_var <- if (length(betas_crude) > 5)
    max(stats::var(betas_crude, na.rm = TRUE), 0.01) else 1.0
  sigma2 <- emp_var * exp(-(0:max_order))
  list(sigma2   = sigma2,
       lambda2  = 1 / pmax(sigma2, 0.01),
       sigma2_g = max(emp_var, 0.25),
       a0 = a0, b0 = b0)
}
