# =============================================================================
# inference.R  --  Posterior inference helpers (cumulative incidence,
#                  interaction effects, etc.)
# =============================================================================


#' Cumulative incidence F_m(tau) for a target node and patient profile
#'
#' Computes \eqn{F_m(\tau) = 1 - \exp(-\int_0^\tau \lambda_m(s)\, ds)}
#' for a given patient trajectory, evaluated under a fitted CTBN. The
#' rate is plugged in at the posterior mean (Stan) or MAP estimate
#' (MAP fit). Works with both \code{ctbn_map_fit} and \code{ctbn_stan_fit}.
#'
#' @param fit A fitted \code{ctbn_fit} object.
#' @param target Character, the target condition.
#' @param patient_dt A data.table of at-risk intervals for one patient.
#' @param tau Numeric, time horizon at which to evaluate F_m.
#' @param zbar Optional named numeric vector to override covariate values.
#' @param unit Optional rate scaling factor.
#'
#' @return Numeric, cumulative incidence in [0, 1].
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' DTw <- prepare_wide(DT)
#' fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 2)
#'
#' # Cumulative incidence of condition "C1" by 5 years for one patient.
#' target  <- fit$call_args$all_conditions[1]
#' pat_dt  <- DTw[eid == DTw$eid[1]]
#' compute_F_m(fit, target = target, patient_dt = pat_dt, tau = 5)
#' }
#'
#' @export
compute_F_m <- function(fit, target, patient_dt, tau,
                          zbar = NULL, unit = 1) {
  ca <- fit$call_args
  all_covs      <- c(ca$fixed_covs, ca$time_varying_covs)
  reserved_base <- c("eid", "time_to_event", "dt", all_covs,
                     paste0(target, "_event"))
  orig_cols     <- names(patient_dt)
  all_conds     <- ca$all_conditions %||% setdiff(orig_cols, reserved_base)
  influencers   <- setdiff(all_conds, target)

  # Get linear predictor on the patient_dt
  pred <- get_lp(fit, patient_dt, target)
  q_hat <- pred$lambda

  if (all(is.na(q_hat))) return(NA_real_)

  # Override covariates if zbar provided -- recompute lp manually.
  if (!is.null(zbar)) {
    nd <- as.data.frame(patient_dt)
    n  <- nrow(nd)
    # Replace covariate columns in patient_dt with zbar values
    for (cv in names(zbar))
      if (cv %in% names(nd)) nd[[cv]] <- zbar[[cv]]
    pred  <- get_lp(fit, data.table::as.data.table(nd), target)
    q_hat <- pred$lambda
  }

  t_start <- as.numeric(patient_dt$time_to_event)
  t_end   <- t_start + as.numeric(patient_dt$dt)
  in_win  <- t_start < tau
  t_star  <- pmax(0, pmin(t_end, tau) - t_start)

  pmin(pmax(1 - exp(-sum(q_hat[in_win] * t_star[in_win], na.rm = TRUE) * unit),
            0), 1)
}


#' Compute pairwise interaction effects (synergistic risk decomposition)
#'
#' For each significant triplet \eqn{(j, k \to m)} where the interaction
#' PIP exceeds \code{pip_thresh}, returns the parameter decomposition
#' (main RRs, synergistic multiplier, joint RR) and the absolute
#' synergistic excess
#' \eqn{\Delta F_m(\tau) = F(j,k) - F(j) - F(k) + F(\emptyset)} at
#' time horizon \code{tau}.
#'
#' @param fit A \code{ctbn_fit} object (MAP or Stan).
#' @param DT_wide Wide-format data.table used for fitting (needed for
#'   covariate reference values).
#' @param tau Numeric, time horizon (years).
#' @param pip_thresh Numeric in [0, 1]; significance threshold for
#'   interaction selection.
#'
#' @return A data.table with columns:
#' \describe{
#'   \item{target_condition, condition_j, condition_k}{names}
#'   \item{beta_j, beta_k, beta_jk}{coefficients}
#'   \item{RR_j, RR_k, syn_multiplier, joint_RR}{rate ratios}
#'   \item{F_baseline, F_j_only, F_k_only, F_joint}{absolute risks}
#'   \item{delta_F}{synergistic excess}
#'   \item{pip_jk}{interaction selection score}
#' }
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500,
#'                            interaction_order = 2, seed = 1)
#' DTw <- prepare_wide(DT)
#'
#' # Interaction decomposition needs a model with max_order >= 2.
#' fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 2)
#' syn <- compute_interaction_effects(fit, DTw, tau = 5, pip_thresh = 0.5)
#' head(syn)
#' }
#'
#' @export
compute_interaction_effects <- function(fit, DT_wide, tau = 5,
                                          pip_thresh = 0.50) {
  ca       <- fit$call_args
  all_covs <- c(ca$fixed_covs, ca$time_varying_covs)
  reserved <- c("eid","time_to_event","dt", all_covs)
  ev_cols  <- grep("_event$", names(DT_wide), value = TRUE)
  all_conds <- setdiff(names(DT_wide), c(reserved, ev_cols))

  if (is.null(ca$max_order) || ca$max_order < 2) {
    method <- ca$method %||% "this method"
    message(sprintf(
      "compute_interaction_effects(): %s fits do not model parent ",
      method),
      "interactions explicitly; returning empty table. Use ",
      "method = 'map' or 'stan' with max_order >= 2 for interaction ",
      "decomposition.")
    return(data.table::data.table())
  }

  # Average covariate profile
  zbar <- vapply(all_covs, function(cv) {
    v <- DT_wide[[cv]]
    if (is.numeric(v)) mean(v, na.rm = TRUE)
    else as.numeric(names(which.max(table(v))))
  }, numeric(1))
  names(zbar) <- all_covs

  results <- list()

  is_map  <- inherits(fit, "ctbn_map_fit")
  is_stan <- inherits(fit, "ctbn_stan_fit")

  for (target in all_conds) {
    # Skip if no model for this target
    has_fit <- if (is_map)  !is.null(fit$map_fits[[target]])
               else if (is_stan) !is.null(fit$stan_fits[[target]])
               else FALSE
    if (!has_fit) next

    influencers <- setdiff(all_conds, target)

    # Locate beta-hat positions for main + 2-way interactions
    if (is_map) {
      mfit <- fit$map_fits[[target]]
      beta_main_mean <- setNames(
        mfit$beta_hat[seq_along(influencers)], influencers)
      x_cols <- mfit$x_cols
      x_orders <- mfit$x_orders
      pos_int2 <- which(x_orders == 1L)
      int_cols <- x_cols[pos_int2]
      beta_int_mean <- setNames(mfit$beta_hat[pos_int2], int_cols)
      pip_int_mean  <- setNames(mfit$pip_hat[pos_int2],  int_cols)
    } else {
      sfit <- fit$stan_fits[[target]]
      post <- as.data.frame(sfit)
      beta_cols <- grep("^beta\\[", names(post), value = TRUE)
      pip_cols  <- grep("^pip_beta\\[", names(post), value = TRUE)

      n_main <- length(influencers)
      tmp <- as.data.frame(matrix(0, 1, n_main))
      colnames(tmp) <- influencers
      ir <- build_interaction_cols(tmp, influencers, ca$max_order)
      int_cols <- ir$cols
      n_int <- length(int_cols)
      if (n_int == 0) next
      if (length(beta_cols) < n_main + n_int) next

      beta_main_mean <- colMeans(post[, beta_cols[seq_len(n_main)],
                                       drop = FALSE])
      names(beta_main_mean) <- influencers
      beta_int_mean <- colMeans(post[, beta_cols[n_main + seq_len(n_int)],
                                      drop = FALSE])
      names(beta_int_mean) <- int_cols
      if (length(pip_cols) >= n_main + n_int) {
        pip_int_mean <- colMeans(post[, pip_cols[n_main + seq_len(n_int)],
                                       drop = FALSE])
      } else {
        pip_int_mean <- setNames(rep(NA_real_, n_int), int_cols)
      }
      names(pip_int_mean) <- int_cols
    }

    for (pair in utils::combn(influencers, 2, simplify = FALSE)) {
      j <- pair[1]; k <- pair[2]
      cn1 <- paste(j, k, sep = "_x_"); cn2 <- paste(k, j, sep = "_x_")
      ac  <- if (cn1 %in% int_cols) cn1
              else if (cn2 %in% int_cols) cn2 else NA
      if (is.na(ac)) next

      beta_jk <- beta_int_mean[[ac]]
      pip_jk  <- pip_int_mean[[ac]]
      beta_j  <- beta_main_mean[[j]]
      beta_k  <- beta_main_mean[[k]]

      if (!is.na(pip_thresh) &&
          (is.na(pip_jk) || pip_jk < pip_thresh)) next

      make_profile <- function(xj, xk) {
        ref <- data.table::as.data.table(
          as.list(rep(0, length(all_conds))))
        data.table::setnames(ref, all_conds)
        ref[, (j)      := as.numeric(xj)]
        ref[, (k)      := as.numeric(xk)]
        ref[, (target) := 0]
        ref[, eid           := 1L]
        ref[, time_to_event := 0]
        ref[, dt            := tau]
        for (cv in all_covs) ref[, (cv) := zbar[[cv]]]
        ref
      }

      F_empty <- tryCatch(compute_F_m(fit, target, make_profile(0,0), tau),
                          error = function(e) NA_real_)
      F_j     <- tryCatch(compute_F_m(fit, target, make_profile(1,0), tau),
                          error = function(e) NA_real_)
      F_k     <- tryCatch(compute_F_m(fit, target, make_profile(0,1), tau),
                          error = function(e) NA_real_)
      F_jk    <- tryCatch(compute_F_m(fit, target, make_profile(1,1), tau),
                          error = function(e) NA_real_)
      delta_F <- F_jk - F_j - F_k + F_empty

      results[[length(results) + 1]] <- data.table::data.table(
        target_condition = target,
        condition_j      = j,
        condition_k      = k,
        beta_j           = round(beta_j,  4),
        beta_k           = round(beta_k,  4),
        beta_jk          = round(beta_jk, 4),
        RR_j             = round(exp(beta_j),  3),
        RR_k             = round(exp(beta_k),  3),
        syn_multiplier   = round(exp(beta_jk), 3),
        joint_RR         = round(exp(beta_j + beta_k + beta_jk), 3),
        F_baseline       = round(F_empty, 4),
        F_j_only         = round(F_j,     4),
        F_k_only         = round(F_k,     4),
        F_joint          = round(F_jk,    4),
        delta_F          = round(delta_F, 4),
        pip_jk           = round(pip_jk,  3),
        tau_years        = tau)
    }
  }

  if (length(results) == 0) return(data.table::data.table())
  data.table::rbindlist(results)
}
