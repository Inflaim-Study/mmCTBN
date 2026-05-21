# =============================================================================
# fit_cph.R  --  CTBN-CPH (Guillamet et al., 2025)
# =============================================================================
#
# Two-stage hybrid (paper, Eq. 11):
#
#   lambda(C_i; m | u) = lambda^0_{m | u}  *  exp( sum_l beta_l * c_l * sigma_l )
#
# Stage 1. Baseline intensities lambda^0_{m|u} come from a *classical*
#          CTBN with the topology supplied (or inferred elsewhere) -- we
#          reuse the closed-form MLE from ctbn_classic().
# Stage 2. For each directed edge (j -> m) in the topology, fit a
#          Cox proportional-hazards model on the trajectory-induced
#          time-to-event for transitions of m, treating patients who
#          have parent j in the relevant state as the at-risk pool.
#          Covariates whose Cox p-value is below `alpha` are kept;
#          others are gated out (sigma_l = 0) when scoring new
#          individuals.
#
# At inference time, lambda for a given (target, parent state, patient
# covariates) is the product of the baseline lambda^0 and the
# patient-specific hazard ratio. This delivers individualised
# trajectories without re-running the structural learner per patient.
# =============================================================================


#' Fit a CTBN-CPH (Guillamet 2025) — individualised CTBN via Cox-PH
#'
#' Hybrid two-stage model. First, baseline rates
#' \eqn{\lambda^0_{m\,|\,u}} are estimated with a classical CTBN (see
#' \code{\link{ctbn_classic}}). Second, for every directed edge
#' \eqn{j \to m} in the topology a Cox proportional-hazards regression
#' is run on the time-to-transition of \eqn{m}, conditioning on
#' patient covariates. At inference time the intensity for a patient
#' with covariates \eqn{C_i} is
#' \deqn{\lambda(C_i;\, m\,|\,u) \;=\; \lambda^0_{m\,|\,u}
#'        \cdot \exp\!\Big( \sum_{l} \beta_l\, c_l\, \sigma_l \Big),}
#' where \eqn{\sigma_l = 1} iff the Cox p-value for covariate \eqn{l}
#' on edge \eqn{j \to m} is below \code{alpha}.
#'
#' @param DT_wide A wide-format data.table; same layout as
#'   \code{\link{ctbn_map}}.
#' @param parents Named list of parent vectors per condition (the
#'   network topology). If \code{NULL}, every other condition is
#'   treated as a parent.
#' @param max_parents Integer cap on per-target parent set size.
#' @param fixed_covs,time_varying_covs Covariate columns to enter the
#'   Cox stage.
#' @param target_conditions Optional subset of conditions to model.
#' @param alpha Numeric significance threshold on Cox p-values; any
#'   covariate with \eqn{p \ge \alpha} is gated out at inference.
#' @param add_laplace Pseudo-count for the classical baseline stage.
#' @param ties Character; ties option passed to
#'   \code{survival::coxph}. Default \code{"efron"}.
#' @param parallel,verbose As elsewhere.
#'
#' @return An object of class \code{c("ctbn_cph_fit", "ctbn_fit")}
#'   with the standard slots plus
#'   \describe{
#'     \item{classic}{The underlying classical CTBN fit (baseline
#'       lambdas).}
#'     \item{cox_fits}{Nested list: \code{cox_fits[[target]][[parent]]}
#'       is a Cox model fitted to time-to-target conditional on parent
#'       presence.}
#'     \item{beta_cov}{Nested list of named numeric vectors of Cox
#'       coefficients per edge.}
#'     \item{p_cov}{Nested list of named numeric vectors of Cox
#'       p-values per edge.}
#'   }
#'
#' @references
#' Guillamet GH, Lopez Segui F, Vidal-Alaball J, Lopez B.
#' \emph{CTBN-PH: A continuous-time Bayesian network for
#' individualised diagnostic risk prediction.} Computers in Biology
#' and Medicine, 197: 111069, 2025.
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' DTw <- prepare_wide(DT)
#'
#' # CTBN-CPH: classical baseline modulated by per-edge Cox regressions.
#' # Requires the 'survival' package.
#' fit <- ctbn_cph(DTw, fixed_covs = c("age", "sex_male"),
#'                  max_parents = 3, alpha = 0.05)
#' summary(fit)
#' }
#'
#' @export
ctbn_cph <- function(DT_wide,
                     parents           = NULL,
                     max_parents       = 4L,
                     fixed_covs        = character(0),
                     time_varying_covs = character(0),
                     target_conditions = NULL,
                     alpha             = 0.05,
                     add_laplace       = 0.5,
                     ties              = "efron",
                     parallel          = FALSE,
                     verbose           = TRUE) {

  if (!requireNamespace("survival", quietly = TRUE))
    stop("ctbn_cph() requires the 'survival' package. ",
         "Install it via install.packages('survival').", call. = FALSE)

  stopifnot(data.table::is.data.table(DT_wide))

  # Stage 1: classical CTBN baseline ------------------------------------------
  if (verbose) message("[ctbn_cph] Stage 1: fitting classical baseline...")
  base <- ctbn_classic(
    DT_wide           = DT_wide,
    parents           = parents,
    max_parents       = max_parents,
    fixed_covs        = fixed_covs,
    time_varying_covs = time_varying_covs,
    target_conditions = target_conditions,
    add_laplace       = add_laplace,
    parallel          = parallel,
    verbose           = verbose)

  conds_all <- base$call_args$all_conditions
  targets   <- base$call_args$target_conditions
  parents   <- base$parents
  all_covs  <- c(fixed_covs, time_varying_covs)

  # Stage 2: per-edge Cox-PH --------------------------------------------------
  if (verbose) message("[ctbn_cph] Stage 2: fitting per-edge Cox-PH models...")

  # Patient-level survival data is built once and reused across edges
  surv_data <- .build_surv_dataset(DT_wide, conds_all = conds_all,
                                    all_covs = all_covs)

  edge_worker <- function(m) {
    res <- list()
    for (j in parents[[m]]) {
      res[[j]] <- .fit_edge_cox(
        surv_data  = surv_data,
        target     = m,
        parent     = j,
        covariates = all_covs,
        ties       = ties,
        verbose    = verbose)
    }
    res
  }

  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    cox_fits <- future.apply::future_lapply(
      targets, edge_worker, future.seed = TRUE)
  } else {
    cox_fits <- lapply(targets, edge_worker)
  }
  names(cox_fits) <- targets

  beta_cov <- lapply(cox_fits, function(per_target) {
    lapply(per_target, function(ef) ef$beta)
  })
  p_cov <- lapply(cox_fits, function(per_target) {
    lapply(per_target, function(ef) ef$pvalue)
  })

  call_args <- c(base$call_args, list(
    method = "cph",
    alpha  = alpha,
    ties   = ties
  ))

  out <- list(
    beta_matrix      = base$beta_matrix,
    pip_matrix       = base$pip_matrix,
    kappa_matrix     = NULL,
    intensity_matrix = base$intensity_matrix,
    classic          = base,
    cox_fits         = cox_fits,
    beta_cov         = beta_cov,
    p_cov            = p_cov,
    parents          = parents,
    call_args        = call_args,
    prior            = "cph"
  )
  class(out) <- c("ctbn_cph_fit", "ctbn_fit")
  out
}


# =============================================================================
# Build a per-patient time-to-event matrix used by the Cox stage
# =============================================================================
# For each (target m, parent j) edge we want to know:
#   - at the moment patient i first has parent j = 1 (entry into the
#     at-risk pool), how long until m turns 1?
# We pre-compute, per patient, the FIRST time each condition flips to 1
# (or NA if never). Then for each (m, j) we form (time, status, covs)
# rows: time = t_first(m) - t_first(j), status = 1 if m occurs after j,
# censored otherwise.
# =============================================================================
.build_surv_dataset <- function(DT, conds_all, all_covs) {

  DT <- data.table::copy(DT)
  data.table::setorder(DT, eid, time_to_event)

  # First-occurrence time per (patient, condition)
  first_time <- list()
  for (cn in conds_all) {
    first_time[[cn]] <- DT[get(cn) == 1, .(t = min(time_to_event)),
                            by = eid]
  }

  # Last-observed time per patient (for censoring)
  last_time <- DT[, .(t_last = max(time_to_event + dt)), by = eid]

  # Per-patient covariate snapshot (mean over trajectory)
  if (length(all_covs)) {
    cov_dt <- DT[, lapply(.SD, function(x) mean(.to_num(x) %||% NA, na.rm = TRUE)),
                  .SDcols = all_covs, by = eid]
  } else {
    cov_dt <- DT[, .(eid = unique(eid))]
  }

  list(first_time = first_time,
       last_time  = last_time,
       cov_dt     = cov_dt,
       conds_all  = conds_all,
       all_covs   = all_covs)
}


# =============================================================================
# Fit one Cox model: time to target m, gated by parent j first-occurrence.
# Patients enter the at-risk pool at t_j, observed until t_m or censoring.
# =============================================================================
.fit_edge_cox <- function(surv_data, target, parent,
                          covariates, ties = "efron",
                          verbose = FALSE) {

  empty <- function() list(
    beta = setNames(rep(0, length(covariates)), covariates),
    pvalue = setNames(rep(1, length(covariates)), covariates),
    n_events = 0L, n_at_risk = 0L, status = "no_data")

  t_par_src <- surv_data$first_time[[parent]]
  if (is.null(t_par_src) || !nrow(t_par_src)) return(empty())
  t_par <- data.table::copy(t_par_src)
  data.table::setnames(t_par, "t", "t_parent")

  t_tgt_src <- surv_data$first_time[[target]]
  if (!is.null(t_tgt_src) && nrow(t_tgt_src)) {
    t_tgt2 <- data.table::copy(t_tgt_src)
    data.table::setnames(t_tgt2, "t", "t_target")
  } else {
    t_tgt2 <- data.table::data.table(eid = integer(0), t_target = numeric(0))
  }

  joined <- merge(t_par, t_tgt2, by = "eid", all.x = TRUE)
  joined <- merge(joined, surv_data$last_time, by = "eid", all.x = TRUE)

  joined[, status := as.integer(!is.na(t_target) & t_target > t_parent)]
  joined[, time := data.table::fifelse(
    status == 1L,
    t_target - t_parent,
    t_last   - t_parent)]
  joined <- joined[time > 0]
  if (!nrow(joined)) return(empty())

  if (length(covariates)) {
    joined <- merge(joined, surv_data$cov_dt, by = "eid", all.x = TRUE)
    # Standardise covariates so the Cox HRs are on a comparable scale
    for (cv in covariates) {
      v <- as.numeric(joined[[cv]])
      v[!is.finite(v)] <- mean(v, na.rm = TRUE)
      sdv <- sd(v, na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      joined[[cv]] <- (v - mean(v, na.rm = TRUE)) / sdv
    }
  }

  # Single-event check
  if (sum(joined$status) < 2L) return(empty())

  fml <- if (length(covariates))
    stats::as.formula(paste(
      "survival::Surv(time, status) ~",
      paste(covariates, collapse = " + ")))
  else
    stats::as.formula("survival::Surv(time, status) ~ 1")

  fit <- tryCatch(
    survival::coxph(fml, data = joined, ties = ties,
                    control = survival::coxph.control(iter.max = 50)),
    error   = function(e) NULL,
    warning = function(w) NULL)

  if (is.null(fit) || length(covariates) == 0L) {
    return(list(
      beta = setNames(rep(0, length(covariates)), covariates),
      pvalue = setNames(rep(1, length(covariates)), covariates),
      n_events = sum(joined$status),
      n_at_risk = nrow(joined),
      status = if (is.null(fit)) "cox_failed" else "no_covariates"))
  }

  smry <- summary(fit)
  beta <- setNames(rep(0, length(covariates)), covariates)
  pval <- setNames(rep(1, length(covariates)), covariates)
  coef_tab <- smry$coefficients
  in_fit   <- intersect(covariates, rownames(coef_tab))
  beta[in_fit] <- coef_tab[in_fit, "coef"]
  pval[in_fit] <- coef_tab[in_fit, "Pr(>|z|)"]
  list(beta = beta, pvalue = pval,
       n_events = sum(joined$status), n_at_risk = nrow(joined),
       status = "ok")
}


# =============================================================================
# get_lp method for CTBN-CPH
# =============================================================================
# For each row of newdata:
#   1. Read the current parent configuration.
#   2. Look up the classical baseline lambda^0 for that configuration.
#   3. For every parent currently in state 1, multiply lambda^0 by
#      exp( sum_l beta^{(j -> m)}_l * c_l * sigma_l ).
#   4. Rows where the target is already present are returned as NA.
# =============================================================================
#' @rdname get_lp
#' @export
get_lp.ctbn_cph_fit <- function(fit, newdata, target, ...) {

  n      <- nrow(newdata)
  base   <- get_lp.ctbn_classic_fit(fit$classic, newdata, target, ...)
  lambda <- base$lambda

  alpha   <- fit$call_args$alpha %||% 0.05
  parents <- fit$parents[[target]]
  if (is.null(parents) || !length(parents)) {
    return(list(lp = log(pmax(lambda, 1e-15)), lambda = lambda))
  }
  covs <- c(fit$call_args$fixed_covs, fit$call_args$time_varying_covs)
  if (!length(covs)) {
    return(list(lp = log(pmax(lambda, 1e-15)), lambda = lambda))
  }

  # Pre-extract covariate matrix from newdata and standardise per
  # patient (matches the standardisation used during Cox fitting).
  C <- matrix(0.0, n, length(covs))
  for (k in seq_along(covs)) {
    v <- as.numeric(.to_num(newdata[[covs[k]]]) %||% rep(0, n))
    v[!is.finite(v)] <- mean(v, na.rm = TRUE)
    sdv <- sd(v, na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    C[, k] <- (v - mean(v, na.rm = TRUE)) / sdv
  }
  colnames(C) <- covs

  hr_log <- numeric(n)
  for (j in parents) {
    parent_state <- as.integer(.to_num(newdata[[j]]) %||% rep(0, n))
    active_rows <- which(parent_state == 1L)
    if (!length(active_rows)) next

    bvec <- fit$beta_cov[[target]][[j]]
    pvec <- fit$p_cov   [[target]][[j]]
    if (is.null(bvec)) next
    sigma <- as.numeric(pvec[covs] < alpha)
    sigma[is.na(sigma)] <- 0
    beta  <- as.numeric(bvec[covs]); beta[is.na(beta)] <- 0
    contrib <- as.numeric(C[active_rows, , drop = FALSE] %*% (beta * sigma))
    hr_log[active_rows] <- hr_log[active_rows] + contrib
  }

  lambda <- lambda * exp(hr_log)
  lp     <- log(pmax(lambda, 1e-15))
  list(lp = lp, lambda = lambda)
}
