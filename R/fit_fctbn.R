# =============================================================================
# fit_fctbn.R  --  Functional CTBN (Faruqui et al., IEEE Access 2021)
# =============================================================================
#
# log( q[x_i | u] ) = beta_0[m] + sum_r gamma_r[m] * z_r
#                                + sum_{j in pa(m)} (
#                                    beta_main[j -> m]
#                                  + sum_r gamma_int[j, r] * z_r
#                                  ) * X_j(t)
#
# Faruqui's central assumption is that *parent effects are
# multiplicative on q*, i.e.\ additive in log q. So the per-parent
# coefficient block is a vector (beta_main + interaction-with-z), and
# higher-order interactions between parents are not modelled separately
# (cf. our MAP/Stan engines, which DO model k-way parent products).
#
# Structure learning uses adaptive group lasso applied per (parent,
# target) coefficient group. Solver: FISTA with a group-soft threshold.
# An optional Gaussian-mixture-based early-stop pushes near-zero groups
# to exactly zero so we don't pay for the long tail of FISTA iterations.
# =============================================================================


#' Fit a Functional CTBN (Faruqui 2021)
#'
#' Models the conditional intensity \eqn{q_{m|u}} as a Poisson
#' regression on exogenous covariates, with parent effects assumed to
#' be multiplicative (i.e.\ additive on the log-rate scale). Structure
#' learning is achieved via adaptive group lasso applied per (parent,
#' target) coefficient block, solved with FISTA. An optional Gaussian
#' mixture model early-stop snaps near-zero blocks to exactly zero,
#' avoiding the long convergence tail of pure proximal gradient.
#'
#' @param DT_wide A wide-format data.table; same layout as
#'   \code{\link{ctbn_map}}.
#' @param fixed_covs,time_varying_covs Character vectors of covariate
#'   column names.
#' @param target_conditions Optional subset of conditions to fit.
#' @param lambda Numeric scalar; the (non-adaptive) base regularisation
#'   strength used to construct adaptive weights as
#'   \eqn{\lambda_j = \lambda / \lVert \tilde\beta_j\rVert} after a
#'   short un-penalised pilot.
#' @param max_iter Integer; FISTA iteration cap. Default 30 000.
#' @param tol Numeric; relative-change convergence tolerance.
#' @param gmm_stop Logical; if TRUE (default), fit a 2-component GMM to
#'   the estimated block norms after a warm-up and zero out groups
#'   inside the near-zero component.
#' @param gmm_warmup Integer; FISTA iterations before the GMM check.
#' @param adaptive Logical; if TRUE (default), apply adaptive group
#'   weights from a pilot unpenalised fit.
#' @param pilot_iter Integer; iterations for the unpenalised pilot.
#' @param parallel,verbose As in \code{\link{ctbn_map}}.
#'
#' @return An object of class \code{c("ctbn_fctbn_fit", "ctbn_fit")}
#'   with the standard slots \code{beta_matrix}, \code{pip_matrix}
#'   (1 - block-soft-threshold indicator), \code{intensity_matrix},
#'   \code{fctbn_fits} (per-target detail) and \code{call_args}.
#'
#' @references
#' Faruqui SHA, Alaeddini A, Wang J, Jaramillo CA, Pugh MJ.
#' \emph{A Functional Model for Structure Learning and Parameter
#' Estimation in Continuous Time Bayesian Network: An Application in
#' Identifying Patterns of Multiple Chronic Conditions.}
#' IEEE Access, 9: 148076–148089, 2021.
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' DTw <- prepare_wide(DT)
#'
#' # Functional CTBN: Poisson regression with adaptive group lasso and
#' # multiplicative parent effects.
#' fit <- ctbn_fctbn(DTw, fixed_covs = c("age", "sex_male"),
#'                    max_iter = 5000)
#' summary(fit)
#' }
#'
#' @export
ctbn_fctbn <- function(DT_wide,
                       fixed_covs        = character(0),
                       time_varying_covs = character(0),
                       target_conditions = NULL,
                       lambda            = 1e-3,
                       max_iter          = 30000L,
                       tol               = 1e-7,
                       gmm_stop          = TRUE,
                       gmm_warmup        = 3000L,
                       adaptive          = TRUE,
                       pilot_iter        = 500L,
                       parallel          = FALSE,
                       verbose           = TRUE) {

  stopifnot(data.table::is.data.table(DT_wide))

  reserved   <- c("eid", "time_to_event", "dt", "fold",
                  fixed_covs, time_varying_covs)
  conds_all  <- setdiff(names(DT_wide),
                        c(reserved,
                          grep("_event$", names(DT_wide), value = TRUE)))

  targets <- target_conditions %||% conds_all
  bad     <- setdiff(targets, conds_all)
  if (length(bad))
    stop("Unknown target conditions: ", paste(bad, collapse = ", "))

  DT <- data.table::copy(DT_wide)
  data.table::setorder(DT, eid, time_to_event)
  for (m in conds_all) {
    ec <- paste0(m, "_event")
    if (!ec %in% names(DT)) {
      DT[, (ec) := as.numeric(
        get(m) == 0 & data.table::shift(get(m), type = "lead") == 1
      ), by = eid]
      DT[is.na(get(ec)), (ec) := 0L]
    }
  }

  worker <- function(m) {
    .fit_fctbn_one_target(
      DT, target = m, conds_all = conds_all,
      fixed_covs = fixed_covs, time_varying_covs = time_varying_covs,
      lambda = lambda, max_iter = max_iter, tol = tol,
      gmm_stop = gmm_stop, gmm_warmup = gmm_warmup,
      adaptive = adaptive, pilot_iter = pilot_iter,
      verbose = verbose)
  }
  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    fits <- future.apply::future_lapply(
      targets, worker, future.seed = TRUE)
  } else {
    fits <- lapply(targets, worker)
  }
  names(fits) <- targets

  bm <- matrix(0, nrow = length(conds_all), ncol = length(conds_all),
               dimnames = list(conds_all, conds_all))
  pm <- matrix(NA_real_, nrow = length(conds_all), ncol = length(conds_all),
               dimnames = list(conds_all, conds_all))
  im <- matrix(NA_real_, nrow = length(conds_all), ncol = length(conds_all),
               dimnames = list(conds_all, conds_all))

  # Intensity at the reference profile (all covariates at their mean ->
  # zero on the standardised scale; all other parents off; the named
  # parent on). With Z's first column being a constant 1, theta[1] is
  # the per-target log-baseline, and theta[groups[[j]]][1] is parent
  # j's main coefficient. So:
  #
  #     log lambda_ref(j -> m) = theta[1] + theta[groups[[j]]][1]
  #
  # which is genuinely distinct from exp(beta_main[j]) because it
  # incorporates the baseline rate.
  for (m in targets) {
    ff   <- fits[[m]]
    bm_m <- ff$beta_main
    pm_m <- ff$selected
    if (is.null(bm_m)) next
    intercept_log <- if (length(ff$theta) >= 1L) ff$theta[1] else 0
    for (j in names(bm_m)) {
      bm[j, m] <- bm_m[[j]]
      pm[j, m] <- as.numeric(pm_m[[j]])
      idx_j    <- ff$groups[[j]]
      beta_j   <- if (length(idx_j) >= 1L) ff$theta[idx_j[1]] else bm_m[[j]]
      im[j, m] <- exp(intercept_log + beta_j)
    }
  }

  call_args <- list(
    method            = "fctbn",
    target_conditions = targets,
    all_conditions    = conds_all,
    fixed_covs        = fixed_covs,
    time_varying_covs = time_varying_covs,
    lambda            = lambda,
    max_iter          = max_iter,
    tol               = tol,
    gmm_stop          = gmm_stop,
    adaptive          = adaptive,
    max_order         = 1L,
    pip_threshold     = 0.5
  )

  out <- list(
    beta_matrix      = bm,
    pip_matrix       = pm,
    kappa_matrix     = NULL,
    intensity_matrix = im,
    fctbn_fits       = fits,
    call_args        = call_args,
    prior            = "fctbn"
  )
  class(out) <- c("ctbn_fctbn_fit", "ctbn_fit")
  out
}


# =============================================================================
# Per-target FCTBN worker
# =============================================================================
.fit_fctbn_one_target <- function(DT, target, conds_all,
                                  fixed_covs, time_varying_covs,
                                  lambda, max_iter, tol,
                                  gmm_stop, gmm_warmup,
                                  adaptive, pilot_iter, verbose) {

  influencers <- setdiff(conds_all, target)
  all_covs    <- c(fixed_covs, time_varying_covs)
  p_cov       <- length(all_covs)

  ec <- paste0(target, "_event")
  atr <- DT[get(target) == 0 & dt > 0]
  n   <- nrow(atr)
  if (n == 0L) {
    return(list(beta_main   = setNames(rep(0, length(influencers)),
                                        influencers),
                selected    = setNames(rep(FALSE, length(influencers)),
                                        influencers),
                groups      = list(),
                theta       = numeric(0),
                ll_path     = numeric(0)))
  }

  y      <- as.integer(atr[[ec]])
  offset <- log(pmax(atr$dt, 1e-12))

  # Build covariate matrix Z (n x (p_cov + 1)) with leading intercept.
  # Standardise non-binary continuous covariates so that the group-lasso
  # penalty applies equivalently across columns and FISTA gradients
  # don't blow up on raw-scale variables like age.
  Z <- matrix(1.0, n, p_cov + 1L)
  Z_centers <- rep(0, p_cov + 1L); Z_centers[1] <- 0
  Z_scales  <- rep(1, p_cov + 1L)
  if (p_cov > 0) {
    for (k in seq_along(all_covs)) {
      v <- .to_num(atr[[all_covs[k]]])
      if (is.null(v)) v <- rep(0, n)
      v[!is.finite(v)] <- mean(v[is.finite(v)], na.rm = TRUE)
      if (!any(is.finite(v))) v <- rep(0, n)
      # Standardise if not effectively binary
      unq <- unique(round(v, 8))
      if (length(unq) > 2L) {
        mu_v  <- mean(v, na.rm = TRUE)
        sd_v  <- stats::sd(v, na.rm = TRUE)
        if (!is.finite(sd_v) || sd_v < 1e-8) sd_v <- 1
        Z[, k + 1L] <- (v - mu_v) / sd_v
        Z_centers[k + 1L] <- mu_v
        Z_scales [k + 1L] <- sd_v
      } else {
        Z[, k + 1L] <- v
      }
    }
  }

  # Per parent j we form a "covariate-modulated indicator" block:
  #   block_j has columns Z .* X_j  (so the first column is just X_j,
  #   subsequent columns are X_j*z_r).
  # This gives parent-specific covariate effects (the multiplicative
  # parent assumption in Faruqui Eq. 9).
  block_cols <- ncol(Z)
  G <- length(influencers)
  Phi <- matrix(0.0, n, ncol(Z) + G * block_cols)
  Phi[, seq_len(ncol(Z))] <- Z

  group_starts <- integer(G + 1L)
  group_starts[1] <- ncol(Z) + 1L
  for (g in seq_len(G)) {
    cn <- influencers[g]
    x  <- .to_num(atr[[cn]]); if (is.null(x)) x <- rep(0, n)
    cols <- group_starts[g]:(group_starts[g] + block_cols - 1L)
    Phi[, cols] <- Z * x
    group_starts[g + 1L] <- group_starts[g] + block_cols
  }
  groups <- lapply(seq_len(G), function(g)
    group_starts[g]:(group_starts[g + 1L] - 1L))
  names(groups) <- influencers

  # ── FISTA configuration ──
  p_total <- ncol(Phi)
  theta   <- rep(0.0, p_total)

  # Lipschitz bound: max eigenvalue of (1/n) * Phi^T diag(mu_max) Phi
  # In the absence of a tight bound, use ||Phi||_F^2 / n.
  L_const <- max(1e-6, sum(Phi * Phi) / n)

  obj_path <- numeric(0)
  prev_obj <- Inf
  theta_prev <- theta
  t_k <- 1
  z_k <- theta

  # Adaptive weights from a pilot (un-penalised) fit. Robust to a
  # divergent / NaN pilot via fallbacks.
  w_g <- rep(1.0, G)
  if (adaptive) {
    pilot <- tryCatch(
      .fista_run(
        Phi = Phi, y = y, offset = offset,
        groups = groups, weights = rep(0.0, G),
        L_const = L_const, max_iter = pilot_iter, tol = tol,
        verbose = verbose, tag = paste0(target, ":pilot")),
      error = function(e) {
        if (verbose)
          message(sprintf(
            "[fctbn:%s] adaptive pilot failed (%s); using uniform weights.",
            target, conditionMessage(e)))
        rep(0.0, ncol(Phi))
      })
    pilot_norms <- vapply(groups, function(idx)
      sqrt(sum(pilot[idx]^2)), numeric(1))
    pilot_norms[!is.finite(pilot_norms)] <- 0
    w_g <- 1 / pmax(pilot_norms, 1e-6)
    w_g[!is.finite(w_g)] <- 1
    mn <- mean(w_g, na.rm = TRUE)
    if (is.finite(mn) && mn > 0) w_g <- w_g / mn else w_g <- rep(1.0, G)
  }

  # ── Main FISTA loop with optional GMM early-stop ──
  res <- .fista_run(
    Phi = Phi, y = y, offset = offset,
    groups = groups, weights = lambda * w_g,
    L_const = L_const, max_iter = max_iter, tol = tol,
    gmm_check = gmm_stop, gmm_warmup = gmm_warmup,
    verbose = verbose, tag = target)

  theta_hat <- res
  theta_hat[!is.finite(theta_hat)] <- 0

  beta_main <- setNames(numeric(G), influencers)
  selected  <- setNames(logical(G),  influencers)
  for (g in seq_len(G)) {
    blk  <- theta_hat[groups[[g]]]
    nrm  <- sqrt(sum(blk^2))
    # First entry of each block is the parent's main coefficient
    # (since Z's first column is the intercept = 1).
    beta_main[g] <- blk[1]
    selected [g] <- nrm > 1e-8
  }

  list(beta_main = beta_main,
       selected  = selected,
       groups    = groups,
       theta     = theta_hat,
       Z_ncol    = ncol(Z),
       cov_names = all_covs,
       Z_centers = Z_centers,
       Z_scales  = Z_scales)
}


# =============================================================================
# Poisson NLL and gradient -- MEAN form for FISTA's convergence test.
#
# These are deliberately namespaced under '.fctbn_' to avoid colliding
# with .poisson_loglik / .poisson_grad in priors.R, which use the
# SUM-form POSITIVE log-likelihood (and gradient) consumed by ctbn_map.
# Conflating the two is a heisenbug: alphabetical source load order
# means priors.R wins, and FISTA then sees a gradient that is (a) the
# wrong sign and (b) ~n times too large, causing immediate divergence.
# =============================================================================
.fctbn_nll <- function(theta, Phi, y, offset) {
  eta <- as.numeric(Phi %*% theta) + offset
  mu  <- exp(pmin(eta, 30))
  sum(mu - y * eta) / length(y)
}

.fctbn_grad <- function(theta, Phi, y, offset) {
  eta <- as.numeric(Phi %*% theta) + offset
  mu  <- exp(pmin(eta, 30))
  as.numeric(crossprod(Phi, mu - y)) / length(y)
}


# =============================================================================
# Group-soft-thresholding proximal operator
# =============================================================================
.group_soft_threshold <- function(theta, groups, lambdas) {
  out <- theta
  for (g in seq_along(groups)) {
    idx <- groups[[g]]
    blk <- theta[idx]
    nrm <- sqrt(sum(blk^2))
    if (nrm <= lambdas[g]) {
      out[idx] <- 0
    } else {
      out[idx] <- (1 - lambdas[g] / nrm) * blk
    }
  }
  out
}


# =============================================================================
# FISTA solver with optional GMM-based early-stop
# =============================================================================
# =============================================================================
# FISTA solver with backtracking line search + optional GMM-based early-stop
#
# Why backtracking. The Poisson loss has Hessian (1/n)*Phi'*diag(mu)*Phi
# where mu = exp(eta), so the true Lipschitz constant depends on the
# current iterate. A static bound like ||Phi||_F^2/n (which is correct
# for least-squares) is far too optimistic for Poisson and lets FISTA
# diverge: eta -> +/-Inf, mu -> exp(30), obj -> Inf or NaN, and the
# convergence test 'abs(prev - obj) / (abs(prev) + eps) < tol' becomes
# NA. Beck & Teboulle (2009)'s backtracking scheme adapts L on the fly
# and is the industry-standard fix.
# =============================================================================
.fista_run <- function(Phi, y, offset, groups, weights,
                       L_const, max_iter, tol = 1e-7,
                       gmm_check = FALSE, gmm_warmup = 3000L,
                       verbose = FALSE, tag = "",
                       L_increase = 2.0, max_bt_steps = 50L) {

  p <- ncol(Phi)
  theta <- rep(0.0, p)
  z     <- theta
  t_k   <- 1.0
  prev  <- .Machine$double.xmax       # finite sentinel (not Inf)
  L     <- max(L_const, 1e-6)

  reg_norm <- function(th)
    sum(weights * vapply(groups, function(i)
      sqrt(sum(th[i]^2)), numeric(1)))

  f_at_z   <- .fctbn_nll(z, Phi, y, offset)
  if (!is.finite(f_at_z)) f_at_z <- .Machine$double.xmax

  for (it in seq_len(max_iter)) {
    g_z <- .fctbn_grad(z, Phi, y, offset)
    if (!all(is.finite(g_z))) {
      # numerical blow-up; reset and bump L
      L     <- L * L_increase * L_increase
      theta <- rep(0.0, p); z <- theta; t_k <- 1.0
      f_at_z <- .fctbn_nll(z, Phi, y, offset)
      next
    }

    # ── Backtracking: find an L such that the proximal step gives an
    # objective <= the local quadratic upper bound at z.
    bt_ok <- FALSE
    for (bt in seq_len(max_bt_steps)) {
      grad_step <- z - g_z / L
      lams      <- weights / L
      theta_new <- .group_soft_threshold(grad_step, groups, lams)

      f_new <- .fctbn_nll(theta_new, Phi, y, offset)
      # Quadratic upper bound Q_L(theta_new, z)
      d <- theta_new - z
      Q <- f_at_z + sum(g_z * d) + 0.5 * L * sum(d * d)

      if (is.finite(f_new) && f_new <= Q + 1e-10) {
        bt_ok <- TRUE
        break
      }
      L <- L * L_increase
      if (L > 1e18) break
    }
    if (!bt_ok) {
      if (verbose)
        message(sprintf(
          "[fctbn:%s] backtracking failed at iter %d (L=%.3g); stopping.",
          tag, it, L))
      break
    }

    # FISTA momentum update
    t_new <- 0.5 * (1 + sqrt(1 + 4 * t_k * t_k))
    z     <- theta_new + ((t_k - 1) / t_new) * (theta_new - theta)
    theta <- theta_new
    t_k   <- t_new
    f_at_z <- .fctbn_nll(z, Phi, y, offset)
    if (!is.finite(f_at_z)) {
      # Drift in the extrapolated point; fall back to non-momentum step
      z      <- theta
      f_at_z <- .fctbn_nll(z, Phi, y, offset)
      t_k    <- 1.0
    }

    # ── Convergence check (safe against Inf/NaN) ──
    if (it %% 100L == 0L || it == max_iter) {
      obj <- f_new + reg_norm(theta)
      if (!is.finite(obj) || !is.finite(prev)) {
        prev <- obj
      } else {
        rel_change <- abs(prev - obj) / (abs(prev) + 1e-9)
        if (is.finite(rel_change) && rel_change < tol) {
          if (verbose)
            message(sprintf(
              "[fctbn:%s] FISTA converged at iter %d (rel.change %.2e)",
              tag, it, rel_change))
          break
        }
        prev <- obj
      }
    }

    if (gmm_check && it == gmm_warmup) {
      g_norms <- vapply(groups, function(i)
        sqrt(sum(theta[i]^2)), numeric(1))
      mask <- .gmm_zero_mask(g_norms)
      for (g in which(mask)) theta[groups[[g]]] <- 0
      z      <- theta
      f_at_z <- .fctbn_nll(z, Phi, y, offset)
      t_k    <- 1.0
      if (verbose && any(mask))
        message(sprintf(
          "[fctbn:%s] GMM early-stop zeroed %d / %d groups at iter %d",
          tag, sum(mask), length(groups), it))
    }
  }
  theta
}


# Two-component GMM (mean ~ 0 vs mean > 0) on group norms; return mask
# of groups in the near-zero component within +/- 3 sd of mu ~ 0.
.gmm_zero_mask <- function(x, eps = 1e-12) {
  if (length(x) < 4 || sd(x) < eps) return(rep(FALSE, length(x)))
  # Initialise: split at the median
  med <- median(x)
  z   <- as.integer(x > med)
  for (iter in seq_len(50)) {
    mu0 <- mean(x[z == 0]); sd0 <- max(sd(x[z == 0]), eps)
    mu1 <- mean(x[z == 1]); sd1 <- max(sd(x[z == 1]), eps)
    if (!is.finite(mu0) || !is.finite(mu1)) break
    p0 <- dnorm(x, mu0, sd0)
    p1 <- dnorm(x, mu1, sd1)
    z_new <- as.integer(p1 > p0)
    if (identical(z_new, z)) break
    z <- z_new
  }
  mu0 <- mean(x[z == 0]); sd0 <- max(sd(x[z == 0]), eps)
  mu1 <- mean(x[z == 1]); sd1 <- max(sd(x[z == 1]), eps)
  # The "near-zero" cluster is whichever has the smaller mean
  zero_cluster <- if (mu0 < mu1) 0 else 1
  mu_z <- if (zero_cluster == 0) mu0 else mu1
  sd_z <- if (zero_cluster == 0) sd0 else sd1
  abs(x - mu_z) <= 3 * sd_z & z == zero_cluster
}


# =============================================================================
# get_lp method for FCTBN
# =============================================================================
#' @rdname get_lp
#' @export
get_lp.ctbn_fctbn_fit <- function(fit, newdata, target, ...) {
  ff <- fit$fctbn_fits[[target]]
  n  <- nrow(newdata)
  if (is.null(ff)) {
    na <- rep(NA_real_, n)
    return(list(lp = na, lambda = na))
  }

  ca       <- fit$call_args
  all_covs <- ff$cov_names
  # Default standardisation constants if fit predates the standardised
  # build (back-compat with older fits).
  Z_centers <- ff$Z_centers %||% rep(0, length(all_covs) + 1L)
  Z_scales  <- ff$Z_scales  %||% rep(1, length(all_covs) + 1L)

  Z <- matrix(1.0, n, length(all_covs) + 1L)
  if (length(all_covs)) {
    for (k in seq_along(all_covs)) {
      v <- .to_num(newdata[[all_covs[k]]])
      if (is.null(v)) v <- rep(0, n)
      v[!is.finite(v)] <- 0
      Z[, k + 1L] <- (v - Z_centers[k + 1L]) / Z_scales[k + 1L]
    }
  }
  influencers <- setdiff(ca$all_conditions, target)
  Phi <- matrix(0.0, n, ncol(Z) + length(influencers) * ncol(Z))
  Phi[, seq_len(ncol(Z))] <- Z
  pos <- ncol(Z)
  for (cn in influencers) {
    x <- .to_num(newdata[[cn]]); if (is.null(x)) x <- rep(0, n)
    x[!is.finite(x)] <- 0
    cols <- (pos + 1L):(pos + ncol(Z))
    Phi[, cols] <- Z * x
    pos <- pos + ncol(Z)
  }
  if (ncol(Phi) != length(ff$theta)) {
    warning(sprintf(
      "get_lp.ctbn_fctbn_fit [%s]: dim mismatch (%d vs %d).",
      target, ncol(Phi), length(ff$theta)))
    na <- rep(NA_real_, n)
    return(list(lp = na, lambda = na))
  }
  lp <- as.numeric(Phi %*% ff$theta)
  lp[!is.finite(lp)] <- NA_real_
  list(lp = lp, lambda = exp(pmin(lp, 30)))
}
