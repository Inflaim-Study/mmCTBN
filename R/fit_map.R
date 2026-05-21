# =============================================================================
# fit_map.R  --  Fast MAP / Laplace estimator for CTBN models
# =============================================================================
#
# This file implements ctbn_map(): the L-BFGS-B MAP optimiser with Laplace
# posterior covariance. It is mathematically equivalent to the Stan MCMC
# fit at the posterior mode but 50-200x faster.
# =============================================================================


# ── Internal: build MAP objective function (closure) ────────────────────────

.make_map_objective <- function(Phi, y, offset, x_orders,
                                 prior, sigma2, lambda2, sigma2_g,
                                 theta_pen, pi0, spike_var,
                                 hs_tau, hs_c2,
                                 n_beta, n_gamma) {
  bi <- seq_len(n_beta)
  gi <- n_beta + seq_len(n_gamma)

  fn <- function(par) {
    beta  <- par[bi];  gamma <- par[gi]
    ll    <- .poisson_loglik(par, Phi, y, offset)
    lp_b  <- switch(prior,
      structured = .prior_logdens_structured(beta, sigma2,  x_orders, theta_pen),
      lasso      = .prior_logdens_lasso(     beta, lambda2, x_orders, theta_pen),
      spike_slab = .prior_logdens_spikeslab( beta, sigma2,  x_orders, theta_pen,
                                              pi0, spike_var),
      horseshoe  = .prior_logdens_horseshoe( beta, hs_tau, hs_c2, x_orders, theta_pen)
    )
    lp_g  <- -0.5 * sum(gamma^2 / sigma2_g + log(2 * pi * sigma2_g))
    -(ll + lp_b + lp_g)
  }

  gr <- function(par) {
    beta  <- par[bi];  gamma <- par[gi]
    gl    <- .poisson_grad(par, Phi, y, offset)
    gp_b  <- switch(prior,
      structured = .grad_prior_structured(beta, sigma2,  x_orders, theta_pen),
      lasso      = .grad_prior_lasso(     beta, lambda2, x_orders, theta_pen),
      spike_slab = .grad_prior_spikeslab( beta, sigma2,  x_orders, theta_pen,
                                           pi0, spike_var),
      horseshoe  = .grad_prior_horseshoe( beta, hs_tau, hs_c2, x_orders, theta_pen)
    )
    gp_g <- -gamma / sigma2_g
    -(gl + c(gp_b, gp_g))
  }

  list(fn = fn, gr = gr)
}


# ── Internal: Laplace SE computation ─────────────────────────────────────────

.laplace_se <- function(theta_hat, Phi, y, offset, x_orders,
                         prior, sigma2, lambda2, sigma2_g,
                         theta_pen, pi0, spike_var, hs_tau, hs_c2,
                         n_beta, n_gamma) {
  d    <- length(theta_hat)
  bi   <- seq_len(n_beta)
  gi   <- n_beta + seq_len(n_gamma)
  beta <- theta_hat[bi]

  H_lik <- .poisson_neg_hess(theta_hat, Phi, offset)

  hd_beta <- switch(prior,
    structured = -.hess_diag_prior_structured(sigma2,  x_orders, theta_pen),
    lasso      = -.hess_diag_prior_lasso(beta, lambda2, x_orders, theta_pen),
    spike_slab = -.hess_diag_prior_spikeslab(beta, sigma2, x_orders, theta_pen,
                                              pi0, spike_var),
    horseshoe  = -.hess_diag_prior_horseshoe(hs_tau, hs_c2, x_orders, theta_pen)
  )

  pd        <- numeric(d)
  pd[bi]    <- hd_beta
  pd[gi]    <- rep(1.0 / sigma2_g, n_gamma)
  H_full    <- H_lik + diag(pd, d)

  Sig <- tryCatch(
    chol2inv(chol(H_full)),
    error = function(e) chol2inv(chol(H_full + diag(1e-6, d)))
  )
  list(se = sqrt(pmax(diag(Sig), 0)), cov_approx = Sig)
}


# ── Internal: selection statistics ──────────────────────────────────────────

.compute_selection_stats <- function(beta_hat, se_beta, x_orders,
                                      prior, sigma2, theta_pen,
                                      pi0, spike_var, hs_tau, hs_c2, n_obs) {
  kappa <- rep(NA_real_, length(beta_hat))

  pip <- switch(prior,
    spike_slab = {
      pi_j   <- pmin(pmax(pi0 * exp(-theta_pen * x_orders), 1e-10), 1 - 1e-10)
      slab_v <- sigma2[x_orders + 1L]
      phi_s  <- stats::dnorm(beta_hat, 0, sqrt(slab_v))
      phi_k  <- stats::dnorm(beta_hat, 0, sqrt(spike_var))
      (pi_j * phi_s) / (pi_j * phi_s + (1 - pi_j) * phi_k + 1e-300)
    },
    horseshoe = {
      v_hs    <- .hs_eff_var(hs_tau, hs_c2, x_orders, theta_pen)
      kappa[] <- 1.0 / (1.0 + n_obs * v_hs)
      1.0 - kappa
    },
    lasso = {
      # Data-aware shrinkage weight matching the horseshoe semantics:
      # small kappa = signal, large kappa = shrunk. Uses the Wald
      # statistic from the Laplace approximation, so coefficients with
      # |beta/se| >> 0 get kappa near 0 and pip near 1.
      w        <- abs(beta_hat) / pmax(se_beta, 1e-10)
      kappa[]  <- 1.0 / (1.0 + w * w)
      stats::pnorm(w)
    },
    stats::pnorm(abs(beta_hat) / pmax(se_beta, 1e-10))
  )

  list(pip = pip, kappa = kappa)
}

.compute_gamma_selection_stats <- function(gamma_hat, se_gamma, sigma2_g,
                                            prior, pi0, spike_var,
                                            hs_tau, hs_c2, n_obs) {
  n_g         <- length(gamma_hat)
  kappa_gamma <- rep(NA_real_, n_g)

  pip_gamma <- switch(prior,
    spike_slab = {
      pi0_g  <- pmin(pmax(pi0, 1e-10), 1 - 1e-10)
      phi_s  <- stats::dnorm(gamma_hat, 0, sqrt(sigma2_g))
      phi_k  <- stats::dnorm(gamma_hat, 0, sqrt(spike_var))
      (pi0_g * phi_s) / (pi0_g * phi_s + (1 - pi0_g) * phi_k + 1e-300)
    },
    horseshoe = {
      v_g_hs      <- .hs_eff_var(hs_tau, hs_c2, orders = rep(0L, n_g), theta = 0)
      kappa_gamma <- 1.0 / (1.0 + n_obs * v_g_hs)
      1.0 - kappa_gamma
    },
    lasso = {
      w           <- abs(gamma_hat) / pmax(se_gamma, 1e-10)
      kappa_gamma <- 1.0 / (1.0 + w * w)
      stats::pnorm(w)
    },
    stats::pnorm(abs(gamma_hat) / pmax(se_gamma, 1e-10))
  )

  list(pip_gamma = pip_gamma, kappa_gamma = kappa_gamma)
}


# =============================================================================
# Main MAP fitting function
# =============================================================================

#' Fit a CTBN by Maximum-A-Posteriori (MAP) estimation
#'
#' Fast L-BFGS-B optimiser for the joint posterior mode plus a Laplace
#' approximation for posterior standard errors. Provides a 50-200x speed-up
#' over the equivalent Stan MCMC fit (\code{\link{ctbn_stan}}) at the cost
#' of replacing posterior samples with normal approximations centred at
#' the MAP.
#'
#' @param DT_wide A data.table in wide format with columns:
#'   \itemize{
#'     \item \code{eid}            -- patient identifier
#'     \item \code{time_to_event}  -- start of at-risk interval
#'     \item \code{dt}             -- length of at-risk interval
#'     \item one column per condition (0/1 indicator)
#'     \item one column per fixed and time-varying covariate
#'   }
#' @param prior Character; one of \code{"spike_slab"}, \code{"structured"},
#'   \code{"lasso"}, \code{"horseshoe"}.
#' @param max_order Integer in [1, 5]; maximum interaction order
#'   (1 = main effects only, 5 = up to 5-way interactions).
#' @param fixed_covs Character vector of fixed-covariate column names.
#' @param time_varying_covs Character vector of time-varying covariate names.
#' @param target_conditions Optional character subset of conditions to fit
#'   as targets. Default \code{NULL} = fit all.
#' @param variable_select Logical; if \code{TRUE}, gate the intensity matrix
#'   by the selection threshold \code{pip_threshold}.
#' @param pip_threshold Numeric in [0, 1]; default 0.5.
#' @param theta Numeric, order penalty strength (>= 0).
#' @param a0,b0 Hyperparameters of the inverse-gamma slab prior on
#'   \eqn{\sigma^2} (structured / lasso / spike_slab).
#' @param pi0 Base inclusion probability (spike_slab only).
#' @param spike_var Spike variance (spike_slab only).
#' @param tau0 Half-Cauchy global scale (horseshoe only).
#' @param slab_df Slab degrees of freedom (horseshoe only).
#' @param slab_scale Slab scale (horseshoe only).
#' @param sigma2_init,lambda2_init,sigma2_g_init Optional manual overrides
#'   for the empirical-Bayes hyperparameter initialisation.
#' @param lbfgs_maxit Integer, maximum L-BFGS-B iterations.
#' @param compute_se Logical; compute Laplace SEs (default TRUE).
#' @param parallel Logical; if TRUE distribute target nodes across workers
#'   via \pkg{future.apply}. Requires \code{plan(multisession)} or similar
#'   to be set by the caller.
#' @param verbose Logical, print progress messages.
#'
#' @return An object of class \code{c("ctbn_map_fit", "ctbn_fit")} with
#'   slots \code{beta_matrix}, \code{se_matrix}, \code{pip_matrix},
#'   \code{kappa_matrix}, \code{intensity_matrix}, \code{pip_cov_list},
#'   \code{convergence}, \code{map_fits}, \code{call_args}.
#'
#'   The \code{intensity_matrix} entry \code{intensity_matrix[j, m]} is
#'   the per-edge rate evaluated at a single, well-defined reference
#'   profile: parent \eqn{j} is on, every \emph{other} parent of
#'   \eqn{m} is off, and the covariates sit at the training-sample
#'   mean (with the intercept column treated as 1 when no covariates
#'   are present). Formally
#'   \deqn{\lambda^{\,\mathrm{ref}}_{j \to m} \;=\;
#'         \exp\!\bigl(\beta_{j,m}^{\,(\text{main})} \;+\;
#'                     \boldsymbol{\gamma}_m'\bar{\mathbf z}\bigr).}
#'   Interactions are \emph{fully part of the fitted model} but
#'   numerically evaluate to zero here, because at the reference state
#'   every interaction column has the form \eqn{X_j \cdot X_k = 1\cdot
#'   0 = 0}. This matrix is therefore the right object to read off the
#'   isolated marginal effect of each parent; it is \emph{not} the
#'   right object for predicting trajectories at arbitrary states (use
#'   \code{\link{compute_F_m}} or \code{\link{get_lp}} for that, both
#'   of which honour every interaction term).
#'
#' @examples
#' \dontrun{
#' # Simulate from a 6-node network with 2-way interactions
#' net <- make_random_network(n_nodes = 6, max_parents = 2, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500,
#'                            interaction_order = 2, seed = 1)
#' DTw <- prepare_wide(DT)
#'
#' # Fast MAP fit
#' fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 2,
#'                  fixed_covs = c("intercept","sex_male"),
#'                  time_varying_covs = "age")
#' summary(fit)
#' }
#'
#' @seealso \code{\link{ctbn_fit}} for the unified entry point;
#'   \code{\link{ctbn_stan}} for full Bayesian MCMC.
#' @export
ctbn_map <- function(DT_wide,
                      prior             = c("spike_slab", "structured",
                                            "lasso", "horseshoe"),
                      max_order         = 1L,
                      fixed_covs        = character(0),
                      time_varying_covs = character(0),
                      target_conditions = NULL,
                      variable_select   = FALSE,
                      pip_threshold     = 0.5,
                      theta             = 1.0,
                      a0                = 2.0,
                      b0                = 1.0,
                      pi0               = 0.5,
                      spike_var         = 0.01,
                      tau0              = 1.0,
                      slab_df           = 4.0,
                      slab_scale        = 2.0,
                      sigma2_init       = NULL,
                      lambda2_init      = NULL,
                      sigma2_g_init     = NULL,
                      lbfgs_maxit       = 500L,
                      compute_se        = TRUE,
                      parallel          = FALSE,
                      verbose           = TRUE) {

  prior <- match.arg(prior)
  stopifnot(data.table::is.data.table(DT_wide))
  stopifnot("eid"           %in% names(DT_wide))
  stopifnot("time_to_event" %in% names(DT_wide))
  stopifnot("dt"            %in% names(DT_wide))
  if (max_order < 1L || max_order > 5L)
    stop("max_order must be an integer in [1, 5].")

  all_covs   <- c(fixed_covs, time_varying_covs)
  all_conds  <- .find_condition_cols(DT_wide,
                                      fixed_covs = fixed_covs,
                                      time_varying_covs = time_varying_covs)
  n_cond <- length(all_conds)
  if (n_cond < 2) stop("Need at least 2 condition columns.")

  target_loop <- if (is.null(target_conditions)) all_conds else {
    bad <- setdiff(target_conditions, all_conds)
    if (length(bad))
      stop("target_conditions not found: ", paste(bad, collapse = ", "))
    target_conditions
  }

  # Horseshoe effective hyperparameters
  hs_tau <- tau0
  hs_c2  <- .hs_c2_mean(slab_df, slab_scale)

  if (verbose) {
    message(sprintf("ctbn_map [%s prior, max_order=%d]: fitting %d targets",
                    prior, max_order, length(target_loop)))
  }

  # Empirical-Bayes initialisation
  hp <- .init_hyperparams_eb(DT_wide, all_conds, all_covs, max_order, a0, b0)
  sigma2   <- rep_len(if (!is.null(sigma2_init))  sigma2_init  else hp$sigma2,
                      max_order + 1L)
  lambda2  <- rep_len(if (!is.null(lambda2_init)) lambda2_init else hp$lambda2,
                      max_order + 1L)
  sigma2_g <- if (!is.null(sigma2_g_init)) sigma2_g_init else hp$sigma2_g

  # Compute event indicators if absent
  DT <- data.table::copy(DT_wide)
  data.table::setorder(DT, eid, time_to_event)
  for (cond in all_conds) {
    ec <- paste0(cond, "_event")
    if (!ec %in% names(DT)) {
      DT[, (ec) := as.numeric(
        get(cond) == 0 & data.table::shift(get(cond), type = "lead") == 1),
        by = eid]
      DT[is.na(get(ec)), (ec) := 0]
    }
  }

  # Per-target worker
  fit_one_target <- function(target) {
    .fit_map_one_target(
      target            = target,
      DT                = DT,
      all_conds         = all_conds,
      all_covs          = all_covs,
      max_order         = max_order,
      prior             = prior,
      sigma2            = sigma2,
      lambda2           = lambda2,
      sigma2_g          = sigma2_g,
      theta_pen         = theta,
      pi0               = pi0,
      spike_var         = spike_var,
      hs_tau            = hs_tau,
      hs_c2             = hs_c2,
      lbfgs_maxit       = lbfgs_maxit,
      compute_se        = compute_se,
      verbose           = verbose
    )
  }

  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE))
      stop("parallel=TRUE requires the 'future.apply' package.")
    fit_results <- future.apply::future_lapply(
      target_loop, fit_one_target,
      future.seed = TRUE,
      future.packages = "data.table")
  } else {
    fit_results <- lapply(target_loop, fit_one_target)
  }
  names(fit_results) <- target_loop

  # Assemble matrices
  mk_mat <- function(fill) matrix(fill, n_cond, n_cond,
                                   dimnames = list(all_conds, all_conds))
  beta_matrix      <- mk_mat(0)
  se_matrix        <- mk_mat(NA_real_)
  pip_matrix       <- mk_mat(NA_real_)
  kappa_matrix     <- mk_mat(NA_real_)
  intensity_matrix <- mk_mat(0)
  convergence_vec  <- setNames(rep(NA_integer_, n_cond), all_conds)
  map_fits         <- setNames(vector("list", n_cond), all_conds)
  pip_cov_list     <- setNames(vector("list", n_cond), all_conds)

  for (target in target_loop) {
    res <- fit_results[[target]]
    if (is.null(res)) next
    influencers <- setdiff(all_conds, target)

    bh <- res$beta_hat; seh <- res$se_hat
    ph <- res$pip_hat;  kh  <- res$kappa_hat

    for (idx in seq_along(influencers)) {
      inf <- influencers[idx]
      beta_matrix[inf,  target] <- bh[idx]
      se_matrix[inf,    target] <- seh[idx]
      pip_matrix[inf,   target] <- ph[idx]
      kappa_matrix[inf, target] <- kh[idx]
    }

    # Intensity at reference
    n_beta <- res$n_beta
    for (idx in seq_along(influencers)) {
      inf <- influencers[idx]
      if (variable_select) {
        pv <- pip_matrix[inf, target]
        kv <- kappa_matrix[inf, target]
        skip <- switch(prior,
          spike_slab = !is.na(pv) && pv < pip_threshold,
          structured = !is.na(pv) && pv < pip_threshold,
          lasso      = !is.na(kv) && kv > (1 - pip_threshold),
          horseshoe  = !is.na(kv) && kv > (1 - pip_threshold),
          FALSE)
        if (isTRUE(skip)) next
      }
      x_ref   <- rep(0.0, n_beta)
      inf_pos <- which(res$x_cols == inf)
      if (length(inf_pos) == 1L) x_ref[inf_pos] <- 1.0
      z_ref <- res$z_ref
      intensity_matrix[inf, target] <- exp(
        sum(x_ref * bh) + sum(z_ref * res$gamma_hat))
    }

    map_fits[[target]]      <- res
    convergence_vec[target] <- res$convergence
    pip_cov_list[[target]]  <- setNames(res$pip_gamma_hat, res$z_cols)
  }

  structure(
    list(
      method           = "map_lbfgs",
      prior            = prior,
      beta_matrix      = beta_matrix,
      se_matrix        = se_matrix,
      pip_matrix       = pip_matrix,
      kappa_matrix     = kappa_matrix,
      intensity_matrix = intensity_matrix,
      pip_cov_list     = pip_cov_list,
      convergence      = convergence_vec,
      map_fits         = map_fits,
      models           = map_fits,
      stan_fits        = setNames(vector("list", n_cond), all_conds),
      pvalue_matrix    = mk_mat(NA_real_),
      call_args        = list(
        method            = "map",
        prior             = prior,
        max_order         = max_order,
        fixed_covs        = fixed_covs,
        time_varying_covs = time_varying_covs,
        target_conditions = target_conditions,
        all_conditions    = all_conds,
        variable_select   = variable_select,
        pip_threshold     = pip_threshold,
        theta             = theta,
        a0 = a0, b0 = b0,
        pi0 = pi0, spike_var = spike_var,
        tau0 = tau0, slab_df = slab_df, slab_scale = slab_scale,
        hs_tau = hs_tau, hs_c2 = hs_c2,
        sigma2 = sigma2, lambda2 = lambda2, sigma2_g = sigma2_g,
        lbfgs_maxit = lbfgs_maxit
      )
    ),
    class = c("ctbn_map_fit", "ctbn_fit")
  )
}


# ── Per-target MAP worker (non-exported) ──────────────────────────────────────

.fit_map_one_target <- function(target, DT, all_conds, all_covs, max_order,
                                 prior, sigma2, lambda2, sigma2_g, theta_pen,
                                 pi0, spike_var, hs_tau, hs_c2,
                                 lbfgs_maxit, compute_se, verbose) {
  ec  <- paste0(target, "_event")
  atr <- DT[get(target) == 0 & dt > 0]
  if (nrow(atr) == 0) {
    warning(sprintf("No at-risk rows for '%s' -- skipping.", target))
    return(NULL)
  }
  df <- as.data.frame(atr)
  n  <- nrow(df)

  # Z matrix
  if (length(all_covs) > 0) {
    Z_mat <- do.call(cbind, lapply(all_covs, function(cv) .to_num(df[[cv]])))
    colnames(Z_mat) <- all_covs
  } else {
    Z_mat <- matrix(1.0, n, 1, dimnames = list(NULL, "__intercept__"))
  }

  # X matrix
  X_all <- do.call(cbind, lapply(all_conds, function(cn) as.numeric(df[[cn]])))
  colnames(X_all) <- all_conds
  dm       <- build_design_matrix(X_all, Z_mat, all_conds, target, max_order)
  Phi      <- dm$Phi
  x_cols   <- dm$x_cols
  x_orders <- dm$x_orders
  n_beta   <- dm$n_beta
  n_gamma  <- dm$n_gamma

  y      <- as.integer(df[[ec]])
  offset <- log(pmax(df$dt, 1e-10))

  obj <- .make_map_objective(
    Phi = Phi, y = y, offset = offset, x_orders = x_orders,
    prior = prior, sigma2 = sigma2, lambda2 = lambda2,
    sigma2_g = sigma2_g, theta_pen = theta_pen,
    pi0 = pi0, spike_var = spike_var,
    hs_tau = hs_tau, hs_c2 = hs_c2,
    n_beta = n_beta, n_gamma = n_gamma)

  opt <- tryCatch(
    stats::optim(par = rep(0.0, n_beta + n_gamma),
                  fn = obj$fn, gr = obj$gr,
                  method  = "L-BFGS-B",
                  control = list(maxit = lbfgs_maxit, factr = 1e7, pgtol = 1e-5)),
    error = function(e) {
      warning(sprintf("L-BFGS-B failed for '%s': %s", target, e$message))
      list(par = rep(0.0, n_beta + n_gamma), convergence = 9L, value = NA)
    })

  theta_hat <- opt$par
  beta_hat  <- theta_hat[seq_len(n_beta)]
  gamma_hat <- theta_hat[n_beta + seq_len(n_gamma)]

  # Laplace SEs
  se_hat   <- rep(NA_real_, n_beta)
  se_gamma <- rep(NA_real_, n_gamma)
  if (compute_se) {
    lap <- tryCatch(
      .laplace_se(theta_hat = theta_hat, Phi = Phi, y = y, offset = offset,
                   x_orders = x_orders, prior = prior,
                   sigma2 = sigma2, lambda2 = lambda2, sigma2_g = sigma2_g,
                   theta_pen = theta_pen, pi0 = pi0, spike_var = spike_var,
                   hs_tau = hs_tau, hs_c2 = hs_c2,
                   n_beta = n_beta, n_gamma = n_gamma),
      error = function(e) {
        warning(sprintf("Laplace SE failed for '%s': %s", target, e$message))
        list(se = rep(NA_real_, n_beta + n_gamma))
      })
    se_hat   <- lap$se[seq_len(n_beta)]
    se_gamma <- lap$se[n_beta + seq_len(n_gamma)]
  }

  sel <- .compute_selection_stats(
    beta_hat = beta_hat, se_beta = se_hat, x_orders = x_orders,
    prior = prior, sigma2 = sigma2, theta_pen = theta_pen,
    pi0 = pi0, spike_var = spike_var,
    hs_tau = hs_tau, hs_c2 = hs_c2, n_obs = n)

  gsel <- .compute_gamma_selection_stats(
    gamma_hat = gamma_hat, se_gamma = se_gamma,
    sigma2_g  = sigma2_g, prior = prior,
    pi0 = pi0, spike_var = spike_var,
    hs_tau = hs_tau, hs_c2 = hs_c2, n_obs = n)

  z_ref <- if (ncol(Z_mat) == 1L && colnames(Z_mat)[1L] == "__intercept__")
    1.0 else colMeans(Z_mat)

  list(
    theta_hat       = theta_hat,
    beta_hat        = beta_hat,
    gamma_hat       = gamma_hat,
    se_hat          = se_hat,
    se_gamma        = se_gamma,
    pip_hat         = sel$pip,
    kappa_hat       = sel$kappa,
    pip_gamma_hat   = gsel$pip_gamma,
    kappa_gamma_hat = gsel$kappa_gamma,
    x_cols          = x_cols,
    x_orders        = x_orders,
    z_cols          = colnames(Z_mat),
    z_ref           = z_ref,
    convergence     = opt$convergence,
    optim_value     = opt$value,
    n_obs           = n,
    n_beta          = n_beta,
    n_gamma         = n_gamma
  )
}
