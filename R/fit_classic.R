# =============================================================================
# fit_classic.R  --  Classical CTBN (Nodelman 2002)
# =============================================================================
#
# Fully parametric continuous-time Bayesian network with discrete state
# variables. No covariates. For each target node m and each configuration
# u of m's parents, the transition rate from "absent" -> "present" is
# estimated by closed-form MLE:
#
#     lambda_hat[m | u]  =  M[m | u]  /  T[m | u]
#
# where M[m|u] is the number of 0->1 transitions for m while its parents
# are in configuration u, and T[m|u] is the total at-risk time spent in
# state 0 while parents are in u. Same for 1 -> 0 transitions.
#
# This module stores rates in a list-of-arrays structure (one array per
# target) and exposes get_lp.ctbn_classic_fit() so the rest of the
# package (CV, metrics, plotting) works uniformly.
# =============================================================================


#' Fit a classical Continuous-Time Bayesian Network (Nodelman 2002)
#'
#' Fully parametric CTBN with discrete (0/1) state variables and no
#' covariates. For every target node and every joint configuration of
#' its parents the 0->1 and 1->0 transition rates are estimated by
#' closed-form maximum-likelihood:
#' \deqn{\hat{\lambda}_{m\,|\,u} \;=\; M[m\,|\,u] \,/\, T[m\,|\,u].}
#'
#' This is the original Nodelman formulation. It does not condition the
#' rates on exogenous covariates and grows exponentially in the number
#' of parents per node — keep \code{max_parents} small.
#'
#' @param DT_wide A wide-format data.table; same layout as for
#'   \code{\link{ctbn_map}}.
#' @param parents Optional named list: \code{parents[[m]]} is the
#'   character vector of parent condition names for target \code{m}. If
#'   \code{NULL}, every other condition is treated as a parent (the
#'   fully connected CTBN).
#' @param max_parents Integer; cap on \code{|parents[[m]]|}. The
#'   parameter count per target is \eqn{2^{|parents|+1}}, so values
#'   above 6 are usually impractical.
#' @param fixed_covs,time_varying_covs Currently ignored (classical
#'   CTBNs have no covariate channel). Retained in the signature so
#'   that \code{ctbn_cv()} drivers can pass them uniformly across
#'   methods.
#' @param target_conditions Optional subset of conditions to fit as
#'   targets.
#' @param add_laplace Numeric pseudo-count added to both
#'   \eqn{M[\cdot]} and \eqn{T[\cdot]} to stabilise rate estimates for
#'   parent configurations that are rare or unobserved. Default 0.5.
#' @param parallel,verbose As in \code{\link{ctbn_map}}.
#'
#' @return An object of class \code{c("ctbn_classic_fit", "ctbn_fit")}
#'   with elements
#'   \describe{
#'     \item{rates}{Named list \code{rates[[m]]} = numeric matrix with
#'       one row per parent configuration and columns
#'       \code{lambda_01}, \code{lambda_10}.}
#'     \item{parent_index}{Named list mapping parent-configuration
#'       strings (e.g.\ "0_1_0") to row indices.}
#'     \item{parents}{Echo of the \code{parents} argument.}
#'     \item{beta_matrix}{Influencer-by-target log-RR matrix
#'       (constructed by marginalising parent-configuration rates so
#'       that the rest of the package — plot, summary, metrics — keeps
#'       working).}
#'     \item{pip_matrix, kappa_matrix}{Inclusion proxies; set to 1 for
#'       every (parent, target) edge in the supplied graph and NA
#'       elsewhere.}
#'     \item{call_args}{Echo of the call.}
#'   }
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' DTw <- prepare_wide(DT)
#'
#' # Classical CTBN: fully parametric, no covariates, parents selected
#' # up to max_parents per node.
#' fit <- ctbn_classic(DTw, max_parents = 3)
#' summary(fit)
#' coef(fit, "intensity")
#' }
#'
#' @export
ctbn_classic <- function(DT_wide,
                         parents           = NULL,
                         max_parents       = 4L,
                         fixed_covs        = character(0),
                         time_varying_covs = character(0),
                         target_conditions = NULL,
                         add_laplace       = 0.5,
                         parallel          = FALSE,
                         verbose           = TRUE) {

  stopifnot(data.table::is.data.table(DT_wide))
  stopifnot("eid" %in% names(DT_wide),
            "time_to_event" %in% names(DT_wide),
            "dt" %in% names(DT_wide))

  reserved   <- c("eid", "time_to_event", "dt", "fold",
                  fixed_covs, time_varying_covs)
  conds_all  <- setdiff(names(DT_wide),
                        c(reserved,
                          grep("_event$", names(DT_wide), value = TRUE)))

  targets <- target_conditions %||% conds_all
  bad     <- setdiff(targets, conds_all)
  if (length(bad))
    stop("Unknown target conditions: ",
         paste(bad, collapse = ", "))

  if (is.null(parents)) {
    parents <- lapply(conds_all, function(m) setdiff(conds_all, m))
    names(parents) <- conds_all
  }
  for (m in targets) {
    if (is.null(parents[[m]])) parents[[m]] <- character(0)
    if (length(parents[[m]]) > max_parents) {
      if (verbose)
        message(sprintf(
          "[ctbn_classic] '%s' has %d parents; truncating to %d.",
          m, length(parents[[m]]), max_parents))
      parents[[m]] <- parents[[m]][seq_len(max_parents)]
    }
  }

  # Ensure target_event columns exist (computed once for the whole DT).
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
    .fit_classic_one_target(DT, target = m, parents_m = parents[[m]],
                            add_laplace = add_laplace)
  }

  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    fits <- future.apply::future_lapply(
      targets, worker, future.seed = TRUE)
  } else {
    fits <- lapply(targets, worker)
  }
  names(fits) <- targets

  # Build beta-style matrices for compatibility with the rest of the
  # package. Edge weight = log(lambda_01 marginal in influencer state 1
  # / lambda_01 marginal in influencer state 0), which approximates a
  # log rate ratio for a single-parent effect.
  bm <- matrix(0, nrow = length(conds_all), ncol = length(conds_all),
               dimnames = list(conds_all, conds_all))
  pm <- matrix(NA_real_, nrow = length(conds_all), ncol = length(conds_all),
               dimnames = list(conds_all, conds_all))
  im <- matrix(NA_real_, nrow = length(conds_all), ncol = length(conds_all),
               dimnames = list(conds_all, conds_all))

  for (m in targets) {
    pm_set <- parents[[m]]
    for (j in pm_set) {
      bm[j, m] <- .classic_log_rr(fits[[m]], parent = j)
      pm[j, m] <- 1
      im[j, m] <- .classic_intensity_at_parent(fits[[m]], parent = j)
    }
  }

  call_args <- list(
    method            = "classic",
    parents           = parents,
    max_parents       = max_parents,
    target_conditions = targets,
    all_conditions    = conds_all,
    fixed_covs        = fixed_covs,
    time_varying_covs = time_varying_covs,
    add_laplace       = add_laplace,
    max_order         = 1L,            # classical CTBN is main-effect
    pip_threshold     = 0.5
  )

  out <- list(
    rates            = lapply(fits, function(f) f$rates),
    parent_index     = lapply(fits, function(f) f$parent_index),
    parents          = parents,
    beta_matrix      = bm,
    pip_matrix       = pm,
    kappa_matrix     = NULL,
    intensity_matrix = im,
    classic_fits     = fits,
    call_args        = call_args,
    prior            = "classic"
  )
  class(out) <- c("ctbn_classic_fit", "ctbn_fit")
  out
}


# ── Per-target worker ────────────────────────────────────────────────────────
.fit_classic_one_target <- function(DT, target, parents_m,
                                    add_laplace = 0.5) {

  ec   <- paste0(target, "_event")
  pa_n <- length(parents_m)

  if (pa_n == 0L) {
    # Marginal rate
    at_risk_0 <- DT[get(target) == 0 & dt > 0]
    at_risk_1 <- data.table::copy(DT[get(target) == 1 & dt > 0])
    M01 <- sum(at_risk_0[[ec]], na.rm = TRUE)
    T01 <- sum(at_risk_0$dt,    na.rm = TRUE)
    # 1 -> 0 transitions
    at_risk_1[, .ev10 := as.numeric(
      get(target) == 1 & data.table::shift(get(target), type = "lead") == 0
    ), by = eid]
    M10 <- sum(at_risk_1$.ev10, na.rm = TRUE)
    T10 <- sum(at_risk_1$dt,    na.rm = TRUE)

    rates <- matrix(c(
      (M01 + add_laplace) / (T01 + add_laplace),
      (M10 + add_laplace) / (T10 + add_laplace)
    ), nrow = 1, dimnames = list(NULL, c("lambda_01", "lambda_10")))
    return(list(rates        = rates,
                parent_index = list("." = 1L),
                parents      = parents_m,
                target       = target))
  }

  # Build parent-configuration key
  config_dt <- data.table::copy(
    DT[, c("eid", "dt", target, ec, parents_m), with = FALSE])
  config_dt[, .pa_key := do.call(paste, c(
    lapply(parents_m, function(p) as.integer(get(p))), sep = "_"))]

  all_configs <- do.call(expand.grid,
                          c(replicate(pa_n, 0:1, simplify = FALSE),
                            stringsAsFactors = FALSE))
  config_keys <- apply(all_configs, 1, paste, collapse = "_")

  rates       <- matrix(NA_real_, nrow = length(config_keys), ncol = 2,
                         dimnames = list(NULL, c("lambda_01", "lambda_10")))
  parent_index <- setNames(as.list(seq_along(config_keys)), config_keys)

  for (i in seq_along(config_keys)) {
    ck <- config_keys[i]
    block <- config_dt[.pa_key == ck]
    if (!nrow(block)) {
      rates[i, ] <- c(add_laplace, add_laplace) /
                    (2 * add_laplace + 1e-9)
      next
    }
    blk0 <- block[get(target) == 0 & dt > 0]
    blk1 <- data.table::copy(block[get(target) == 1 & dt > 0])
    M01  <- sum(blk0[[ec]], na.rm = TRUE)
    T01  <- sum(blk0$dt,    na.rm = TRUE)

    blk1[, .ev10 := as.numeric(
      get(target) == 1 & data.table::shift(get(target), type = "lead") == 0
    ), by = eid]
    M10 <- sum(blk1$.ev10, na.rm = TRUE)
    T10 <- sum(blk1$dt,    na.rm = TRUE)

    rates[i, "lambda_01"] <- (M01 + add_laplace) / (T01 + add_laplace)
    rates[i, "lambda_10"] <- (M10 + add_laplace) / (T10 + add_laplace)
  }
  rownames(rates) <- config_keys

  list(rates        = rates,
       parent_index = parent_index,
       parents      = parents_m,
       target       = target)
}


# ── Per-influencer intensity at a baseline parent configuration ─────────────
#
# Returns lambda_01 evaluated at the parent configuration where only
# `parent` is on and every other parent in the target's parent set is
# 0. This is the natural analogue of MAP/Stan's "intensity at the
# reference profile" and is genuinely distinct from RR = exp(beta).
#
# When the requested config is unobserved in the data, the rates matrix
# already carries the Laplace-smoothed pseudo-rate, so this never
# returns NA.
.classic_intensity_at_parent <- function(target_fit, parent) {
  pa_set <- target_fit$parents
  if (!parent %in% pa_set) return(NA_real_)
  rates  <- target_fit$rates
  states <- as.integer(pa_set == parent)        # 1 for `parent`, 0 elsewhere
  key    <- paste(states, collapse = "_")
  idx    <- target_fit$parent_index[[key]]
  if (is.null(idx)) return(NA_real_)
  rates[idx, "lambda_01"]
}

.classic_log_rr <- function(target_fit, parent) {
  pa_set <- target_fit$parents
  if (!parent %in% pa_set) return(0)
  rates    <- target_fit$rates
  keys     <- rownames(rates)
  pa_pos   <- match(parent, pa_set)
  splits   <- strsplit(keys, "_", fixed = TRUE)
  pa_vals  <- vapply(splits, function(s) s[pa_pos], character(1))

  mean01_p1 <- mean(rates[pa_vals == "1", "lambda_01"], na.rm = TRUE)
  mean01_p0 <- mean(rates[pa_vals == "0", "lambda_01"], na.rm = TRUE)
  if (!is.finite(mean01_p0) || mean01_p0 <= 0 ||
      !is.finite(mean01_p1) || mean01_p1 <= 0) return(0)
  log(mean01_p1 / mean01_p0)
}


# ── get_lp method ────────────────────────────────────────────────────────────
#
# For a classical CTBN, lambda(t) depends only on the *current* parent
# configuration -- not on covariates. We return:
#   lambda = rates[config(t), "lambda_01"]
#   lp     = log(lambda)
# rows where the patient already has the target are returned as NA
# (no 0 -> 1 transition possible).

#' @rdname get_lp
#' @export
get_lp.ctbn_classic_fit <- function(fit, newdata, target, ...) {
  cf <- fit$classic_fits[[target]]
  n  <- nrow(newdata)
  if (is.null(cf)) {
    na <- rep(NA_real_, n)
    return(list(lp = na, lambda = na))
  }
  pa_set <- cf$parents
  if (length(pa_set) == 0L) {
    lam <- rep(cf$rates[1, "lambda_01"], n)
    lp  <- log(pmax(lam, 1e-15))
    # NA-out rows already in state 1
    is_one <- as.integer(newdata[[target]] == 1L)
    lam[is_one == 1L] <- NA_real_
    lp [is_one == 1L] <- NA_real_
    return(list(lp = lp, lambda = lam))
  }

  key <- do.call(paste, c(
    lapply(pa_set, function(p) as.integer(newdata[[p]])),
    sep = "_"))

  idx <- match(key, rownames(cf$rates))
  lam <- cf$rates[idx, "lambda_01"]
  lam[is.na(lam)] <- exp(fit$beta_matrix[setdiff(rownames(fit$beta_matrix),
                                                  target)[1], target] *
                          0)  # fallback: 1.0
  is_one <- as.integer(newdata[[target]] == 1L)
  lam[is_one == 1L] <- NA_real_

  lp <- log(pmax(lam, 1e-15))
  list(lp = lp, lambda = lam)
}
