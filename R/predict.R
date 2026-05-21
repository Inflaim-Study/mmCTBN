# =============================================================================
# predict.R  --  get_lp() S3 generic and methods
# =============================================================================
#
# get_lp(fit, newdata, target) returns a list with `lp` (linear predictor
# eta = X*beta + Z*gamma) and `lambda = exp(lp)` for the given target node
# evaluated on newdata. Used downstream for predictive metrics and CV.
# =============================================================================


#' Linear predictor and rate for a CTBN fit
#'
#' Generic that returns the linear predictor \eqn{\eta = X\beta + Z\gamma}
#' and rate \eqn{\lambda = \exp(\eta)} for one target node, evaluated on
#' new data. S3 methods dispatch on whether \code{fit} is a MAP fit
#' (\code{ctbn_map_fit}) or a Stan fit (\code{ctbn_stan_fit}).
#'
#' @param fit A fitted \code{ctbn_fit} object.
#' @param newdata A data.table of at-risk intervals.
#' @param target Character, the target condition.
#' @param ... Currently unused.
#'
#' @return A list with elements \code{lp} (linear predictor) and
#'   \code{lambda} (rate).
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' DTw <- prepare_wide(DT)
#' fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 1)
#'
#' target <- fit$call_args$all_conditions[1]
#' pred   <- get_lp(fit, DTw, target)
#' str(pred)
#' }
#'
#' @export
get_lp <- function(fit, newdata, target, ...) UseMethod("get_lp")


#' @rdname get_lp
#' @export
get_lp.ctbn_map_fit <- function(fit, newdata, target, ...) {
  mfit <- fit$map_fits[[target]]
  if (is.null(mfit)) {
    na <- rep(NA_real_, nrow(newdata))
    return(list(lp = na, lambda = na))
  }

  ca        <- fit$call_args
  all_covs  <- c(ca$fixed_covs, ca$time_varying_covs)
  all_conds <- ca$all_conditions
  nd        <- as.data.frame(newdata)
  n         <- nrow(nd)

  # Drop event columns
  nd <- nd[, setdiff(names(nd),
                       grep("_event$", names(nd), value = TRUE)),
           drop = FALSE]

  X_all <- do.call(cbind, lapply(all_conds, function(cn) {
    .to_num(nd[[cn]]) %||% rep(0, n)
  }))
  colnames(X_all) <- all_conds

  Z_mat <- if (length(all_covs) > 0L)
    do.call(cbind, lapply(all_covs, function(cv) {
      .to_num(nd[[cv]]) %||% rep(0, n)
    }))
  else
    matrix(1.0, n, 1L)

  dm <- build_design_matrix(X_all, Z_mat, all_conds, target, ca$max_order)

  if (ncol(dm$Phi) != (mfit$n_beta + mfit$n_gamma)) {
    warning(sprintf("get_lp [%s]: dim mismatch -- %d vs %d",
                    target, ncol(dm$Phi), mfit$n_beta + mfit$n_gamma))
    na <- rep(NA_real_, n)
    return(list(lp = na, lambda = na))
  }

  lp <- as.numeric(dm$Phi %*% mfit$theta_hat)
  list(lp = lp, lambda = exp(lp))
}


#' @rdname get_lp
#' @export
get_lp.ctbn_stan_fit <- function(fit, newdata, target, ...) {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("get_lp.ctbn_stan_fit() requires 'rstan'.")

  newdata  <- data.table::copy(newdata)
  stan_fit <- fit$stan_fits[[target]]
  if (is.null(stan_fit)) {
    na <- rep(NA_real_, nrow(newdata))
    return(list(lp = na, lambda = na))
  }

  post    <- as.data.frame(stan_fit)
  b_samps <- as.matrix(post[, grep("^beta\\[",  names(post)), drop = FALSE])
  g_samps <- as.matrix(post[, grep("^gamma\\[", names(post)), drop = FALSE])
  beta_mean  <- colMeans(b_samps)
  gamma_mean <- colMeans(g_samps)

  ca       <- fit$call_args
  all_covs <- c(ca$fixed_covs, ca$time_varying_covs)
  all_conds_fit <- ca$all_conditions

  event_cols <- grep("_event$", names(newdata), value = TRUE)
  if (length(event_cols))
    newdata[, (event_cols) := NULL]

  influencers <- setdiff(all_conds_fit, target)
  x_cols      <- influencers
  nd          <- as.data.frame(newdata)
  n           <- nrow(nd)

  if (!is.null(ca$max_order) && ca$max_order >= 2) {
    int_res  <- build_interaction_cols(nd, influencers, ca$max_order)
    nd       <- int_res$data
    x_cols   <- c(x_cols, int_res$cols)
  }

  X_new <- matrix(0, nrow = n, ncol = length(x_cols))
  for (j in seq_along(x_cols))
    X_new[, j] <- .to_num(nd[[x_cols[j]]]) %||% rep(0, n)

  if (length(all_covs) > 0) {
    Z_new <- matrix(0, nrow = n, ncol = length(all_covs))
    for (k in seq_along(all_covs))
      Z_new[, k] <- .to_num(nd[[all_covs[k]]]) %||% rep(0, n)
  } else {
    Z_new <- matrix(1, nrow = n, ncol = 1)
  }

  if (length(beta_mean) != ncol(X_new) ||
      length(gamma_mean) != ncol(Z_new)) {
    warning(sprintf("get_lp [%s]: dim mismatch -- beta %d vs X %d; gamma %d vs Z %d",
                    target, length(beta_mean), ncol(X_new),
                    length(gamma_mean), ncol(Z_new)))
    na <- rep(NA_real_, n)
    return(list(lp = na, lambda = na))
  }

  lp <- as.numeric(X_new %*% beta_mean + Z_new %*% gamma_mean)
  list(lp = lp, lambda = exp(lp))
}


#' @rdname get_lp
#' @export
get_lp.ctbn_fit <- function(fit, newdata, target, ...) {
  if (inherits(fit, "ctbn_map_fit"))     return(get_lp.ctbn_map_fit    (fit, newdata, target, ...))
  if (inherits(fit, "ctbn_stan_fit"))    return(get_lp.ctbn_stan_fit   (fit, newdata, target, ...))
  if (inherits(fit, "ctbn_cph_fit"))     return(get_lp.ctbn_cph_fit    (fit, newdata, target, ...))
  if (inherits(fit, "ctbn_fctbn_fit"))   return(get_lp.ctbn_fctbn_fit  (fit, newdata, target, ...))
  if (inherits(fit, "ctbn_classic_fit")) return(get_lp.ctbn_classic_fit(fit, newdata, target, ...))

  # Heuristic fallback when class info is absent but slot names hint
  if (!is.null(fit$stan_fits)    && !is.null(fit$stan_fits   [[target]]))
    return(get_lp.ctbn_stan_fit    (fit, newdata, target, ...))
  if (!is.null(fit$map_fits)     && !is.null(fit$map_fits    [[target]]))
    return(get_lp.ctbn_map_fit     (fit, newdata, target, ...))
  if (!is.null(fit$cox_fits)     && !is.null(fit$cox_fits    [[target]]))
    return(get_lp.ctbn_cph_fit     (fit, newdata, target, ...))
  if (!is.null(fit$fctbn_fits)   && !is.null(fit$fctbn_fits  [[target]]))
    return(get_lp.ctbn_fctbn_fit   (fit, newdata, target, ...))
  if (!is.null(fit$classic_fits) && !is.null(fit$classic_fits[[target]]))
    return(get_lp.ctbn_classic_fit (fit, newdata, target, ...))

  na <- rep(NA_real_, nrow(newdata))
  list(lp = na, lambda = na)
}
