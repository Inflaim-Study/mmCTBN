# =============================================================================
# metrics.R  --  Simulation-study metrics: oracle, recovery, selection, pred
# =============================================================================
#
# These metrics compare a fitted CTBN model to its known simulation truth,
# and are used by ctbn_cv() and the simulation runner.
# =============================================================================


#' Oracle (true-rate) predictive metrics on test data
#'
#' Computes the Poisson log-likelihood and Brier score using the
#' \emph{true} simulation rates -- i.e., the upper-bound performance
#' achievable if the data-generating process were known.
#'
#' @param DT_test A wide-format data.table including the per-target
#'   event indicators (\code{<cond>_event}).
#' @param network The network spec used to generate the data.
#' @param interaction_order Integer in [1, 5]; truth's interaction order.
#'
#' @return A data.table with columns \code{condition}, \code{oracle_pll},
#'   \code{oracle_brier}.
#'
#' @examples
#' \dontrun{
#' net    <- make_random_network(n_nodes = 6, seed = 1)
#' DTtest <- prepare_wide(
#'   simulate_ctbn_data(network = net, n_patients = 300,
#'                       interaction_order = 2, seed = 2))
#'
#' # Upper-bound performance under the known data-generating process.
#' compute_oracle_metrics(DTtest, network = net, interaction_order = 2)
#' }
#'
#' @export
compute_oracle_metrics <- function(DT_test, network,
                                     interaction_order = 3L) {
  conds <- network$conds
  cov_names <- network$cov_names %||%
    c("age","sex_male","smk_current","smk_former")

  DT_test <- .ensure_event_cols(DT_test, conds)

  results <- list()
  for (m in conds) {
    ec  <- paste0(m, "_event")
    atr <- DT_test[get(m) == 0 & dt > 0]
    if (nrow(atr) == 0 || !ec %in% names(atr)) next
    y   <- atr[[ec]]; dt_ <- atr$dt

    q_true <- vapply(seq_len(nrow(atr)), function(ri) {
      row  <- atr[ri]
      X_st <- setNames(unlist(row[, ..conds]), conds)
      covs <- vapply(cov_names, function(cv) {
        if (cv %in% names(row)) as.numeric(row[[cv]]) else 0
      }, numeric(1))
      names(covs) <- cov_names
      .compute_q_true(network, m, X_st, covs, max_order = interaction_order)
    }, numeric(1))

    mu  <- q_true * pmax(dt_, 1e-10)
    p   <- 1 - exp(-mu)
    results[[m]] <- data.table::data.table(
      condition    = m,
      oracle_pll   = mean(stats::dpois(y, pmax(mu, 1e-15), log = TRUE),
                           na.rm = TRUE),
      oracle_brier = mean((y - p)^2, na.rm = TRUE))
  }
  data.table::rbindlist(results)
}


#' Predictive metrics on test data (Poisson log-likelihood, Brier, AUC)
#'
#' @param fit A fitted \code{ctbn_fit}.
#' @param DT_test Test data.
#' @param eval_times Numeric vector of evaluation horizons for
#'   time-dependent AUC (default \code{1L}).
#'
#' @return A data.table with columns \code{condition}, \code{pll},
#'   \code{brier}, plus one column per \code{eval_times}.
#'
#' @examples
#' \dontrun{
#' net    <- make_random_network(n_nodes = 6, seed = 1)
#' DT     <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' DTw    <- prepare_wide(DT)
#' DTtest <- prepare_wide(
#'   simulate_ctbn_data(network = net, n_patients = 300, seed = 2))
#'
#' fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 1)
#' compute_pred_metrics(fit, DTtest, eval_times = c(2, 5))
#' }
#'
#' @export
compute_pred_metrics <- function(fit, DT_test, eval_times = 1L) {
  ca       <- fit$call_args
  conds    <- ca$all_conditions
  DT_test  <- .ensure_event_cols(DT_test, conds)
  results  <- list()

  for (m in conds) {
    ec  <- paste0(m, "_event")
    atr <- DT_test[get(m) == 0 & dt > 0]
    if (nrow(atr) == 0 || !ec %in% names(atr)) next

    pred <- tryCatch(get_lp(fit, atr, m),
                      error = function(e)
                        list(lp = rep(NA, nrow(atr)),
                             lambda = rep(NA, nrow(atr))))
    y   <- atr[[ec]]; dt_ <- atr$dt
    lam <- pred$lambda; lp <- pred$lp
    mu  <- lam * pmax(dt_, 1e-10)
    p_h <- 1 - exp(-pmax(mu, 0))
    pll <- mean(stats::dpois(y, pmax(mu, 1e-15), log = TRUE), na.rm = TRUE)
    brier <- mean((y - p_h)^2, na.rm = TRUE)

    # ── Per-patient survival table for time-dependent AUC ─────────────
    # Row-level data is the wrong unit of analysis for timeROC: each
    # patient appears multiple times with different intervals, and the
    # at-risk set at large horizons is selectively "patients who
    # haven't transitioned yet" -- which inflates AUC with t. The
    # correct construction is:
    #   T_i      = patient i's first transition time, or last observed
    #              time if censored;
    #   delta_i  = 1 if patient transitioned, 0 otherwise;
    #   marker_i = the patient's baseline log-rate (lp at their first
    #              at-risk row), a prognostic score that does not
    #              depend on the response.
    av <- setNames(rep(NA_real_, length(eval_times)),
                    paste0("auc_t", eval_times))

    if (requireNamespace("timeROC", quietly = TRUE) &&
        "eid" %in% names(atr)) {
      atr_use <- atr[is.finite(lp)]
      lp_use  <- lp[is.finite(lp)]
      if (nrow(atr_use) > 1L) {
        atr_use <- data.table::copy(atr_use)
        atr_use[, `:=`(.row_lp = lp_use,
                        .row_y  = atr_use[[ec]])]
        data.table::setorder(atr_use, eid, time_to_event)

        first_row <- atr_use[, .SD[1L], by = eid]

        # Per-patient event time / censoring
        ev_rows <- atr_use[get(ec) == 1L,
                            .(t_ev = min(time_to_event + dt)),
                            by = eid]
        last_rows <- atr_use[, .(t_last = max(time_to_event + dt)),
                              by = eid]

        pd <- merge(first_row[, .(eid, marker = .row_lp)],
                     last_rows, by = "eid", all.x = TRUE)
        pd <- merge(pd, ev_rows, by = "eid", all.x = TRUE)
        pd[, `:=`(event  = as.integer(!is.na(t_ev)),
                  T_obs  = data.table::fifelse(!is.na(t_ev), t_ev, t_last))]
        pd <- pd[is.finite(marker) & is.finite(T_obs) & T_obs > 0]

        if (nrow(pd) > 1L && length(unique(pd$event)) == 2L) {
          td <- tryCatch(
            timeROC::timeROC(T = pd$T_obs, delta = pd$event,
                              marker = pd$marker,
                              cause = 1, times = eval_times,
                              iid = FALSE),
            error = function(e) NULL)
          if (!is.null(td)) {
            # timeROC sometimes pre-pends t=0 (when length(times)==1)
            # and labels output with "t=<val>" names. Realign by name
            # so we always return exactly one column per requested
            # eval_time in the requested order.
            auc_named <- as.numeric(td$AUC)
            names(auc_named) <- names(td$AUC)
            wanted <- paste0("t=", eval_times)
            picked <- auc_named[wanted]
            av     <- setNames(unname(picked),
                                paste0("auc_t", eval_times))
          }
        }
      }
    }

    results[[m]] <- do.call(data.table::data.table,
      c(list(condition = m, pll = pll, brier = brier), as.list(av)))
  }

  data.table::rbindlist(results, fill = TRUE)
}


#' Selection metrics (TPR / FPR / AUC for variable selection)
#'
#' @param fit A fitted \code{ctbn_fit} with selection scores.
#' @param network The true network specification.
#' @param pip_thresh Selection threshold (default 0.5).
#' @param max_order_truth Highest-order interactions to evaluate
#'   (>= fit's max_order is fine).
#'
#' @return A data.table with one row per interaction order
#'   (\code{main_effect}, \code{interaction_2way}, \code{interaction_3way},
#'   ...) and columns \code{tpr}, \code{fpr}, \code{selection_auc}.
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500,
#'                            interaction_order = 2, seed = 1)
#' DTw <- prepare_wide(DT)
#'
#' fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 2,
#'                  variable_select = TRUE)
#' compute_selection_metrics(fit, network = net, pip_thresh = 0.5)
#' }
#'
#' @export
compute_selection_metrics <- function(fit, network,
                                        pip_thresh = 0.5,
                                        max_order_truth = NULL) {

  if (!inherits(fit, c("ctbn_map_fit", "ctbn_stan_fit"))) {
    method <- fit$call_args$method %||% class(fit)[1]
    message(sprintf(
      "compute_selection_metrics(): not implemented for method '%s' ",
      method),
      "(applies to 'map' and 'stan' only -- those expose per-coefficient ",
      "PIP / kappa scores). Use coef(fit, 'pip') or coef(fit, 'beta') for ",
      "a simpler edge-level overview.")
    return(data.table::data.table())
  }

  ca <- fit$call_args
  if (is.null(max_order_truth)) max_order_truth <- ca$max_order
  is_map <- inherits(fit, "ctbn_map_fit")

  conds <- network$conds
  prior <- fit$prior

  true_main_nz <- names(network$beta_main)[network$beta_main != 0]
  true_2way_nz <- if (length(network$beta_int2))
    names(network$beta_int2)[unlist(network$beta_int2) != 0] else character(0)
  true_3way_nz <- if (length(network$beta_int3))
    names(network$beta_int3)[unlist(network$beta_int3) != 0] else character(0)
  true_4way_nz <- if (length(network$beta_int4))
    names(network$beta_int4)[unlist(network$beta_int4) != 0] else character(0)
  true_5way_nz <- if (length(network$beta_int5))
    names(network$beta_int5)[unlist(network$beta_int5) != 0] else character(0)

  mk_score <- function(mfit, pos) {
    if (is.null(mfit) || pos > length(mfit$beta_hat %||% c()))
      return(NA_real_)
    if (prior == "horseshoe") {
      kh <- mfit$kappa_hat
      if (!is.null(kh) && pos <= length(kh)) return(1.0 - kh[pos])
    } else {
      ph <- mfit$pip_hat
      if (!is.null(ph) && pos <= length(ph)) return(ph[pos])
    }
    NA_real_
  }
  is_selected <- function(score) !is.na(score) && score >= pip_thresh

  # Bucketised rows; we accumulate via direct list-element assignment in
  # the surrounding scope (no closure-based assignment, which would be
  # local to the helper and lost on return).
  rows_main <- list(); rows_int2 <- list(); rows_int3 <- list()
  rows_int4 <- list(); rows_int5 <- list()

  for (m in conds) {
    if (is_map) {
      mfit <- fit$map_fits[[m]]
      if (is.null(mfit)) next
    } else {
      sfit <- fit$stan_fits[[m]]
      if (is.null(sfit)) next
      post <- as.data.frame(sfit)
      beta_cols <- grep("^beta\\[",      names(post), value = TRUE)
      pip_cols  <- grep("^pip_beta\\[",  names(post), value = TRUE)
      kap_cols  <- grep("^kappa\\[",     names(post), value = TRUE)
      mfit <- list(
        beta_hat  = colMeans(post[, beta_cols, drop = FALSE]),
        pip_hat   = if (length(pip_cols))
                      colMeans(post[, pip_cols, drop = FALSE])
                    else rep(NA_real_, length(beta_cols)),
        kappa_hat = if (length(kap_cols))
                      colMeans(post[, kap_cols, drop = FALSE])
                    else rep(NA_real_, length(beta_cols)))
      influencers <- setdiff(conds, m)
      tmp <- as.data.frame(matrix(0, 1, length(influencers)))
      colnames(tmp) <- influencers
      ir <- build_interaction_cols(tmp, influencers, ca$max_order)
      mfit$x_cols   <- c(influencers, ir$cols)
      mfit$x_orders <- c(rep(0L, length(influencers)), ir$orders)
    }

    infl     <- setdiff(conds, m)
    x_cols   <- mfit$x_cols
    x_orders <- mfit$x_orders

    # Main effects -----------------------------------------------------
    for (idx in seq_along(infl)) {
      j   <- infl[idx]
      key <- paste(m, j, sep = "|")
      sc  <- mk_score(mfit, idx)
      rows_main[[length(rows_main) + 1L]] <- data.table::data.table(
        target = m, influencer = j,
        true_nz  = key %in% true_main_nz,
        selected = is_selected(sc), score = sc)
    }

    # Higher-order interactions ----------------------------------------
    # We accumulate rows in named scalars (rows_int2, ...) so that
    # mutation in this for-loop body actually persists -- a closure
    # would silently drop the writes.
    for (ord_int in 2L:5L) {
      if (ord_int > ca$max_order || ord_int > max_order_truth) next

      true_keys <- switch(as.character(ord_int),
        "2" = true_2way_nz, "3" = true_3way_nz,
        "4" = true_4way_nz, "5" = true_5way_nz, character(0))

      grps <- utils::combn(infl, ord_int, simplify = FALSE)
      for (grp in grps) {
        # Try the canonical and all permuted column names; design.R
        # uses combn() order, which is always sorted by appearance in
        # 'influencers', so the first match below should be sufficient.
        perms <- vapply(.permutations(grp),
                          function(p) paste(p, collapse = "_x_"),
                          character(1))
        ac <- intersect(perms, x_cols)
        if (!length(ac)) next
        pos <- which(x_cols == ac[1])[1]
        sc  <- mk_score(mfit, pos)
        keys <- vapply(.permutations(grp),
                        function(p) paste(c(m, p), collapse = "|"),
                        character(1))
        rec <- data.table::data.table(
          target   = m,
          group    = paste(grp, collapse = "_x_"),
          true_nz  = any(keys %in% true_keys),
          selected = is_selected(sc),
          score    = sc)

        if      (ord_int == 2L) rows_int2[[length(rows_int2) + 1L]] <- rec
        else if (ord_int == 3L) rows_int3[[length(rows_int3) + 1L]] <- rec
        else if (ord_int == 4L) rows_int4[[length(rows_int4) + 1L]] <- rec
        else if (ord_int == 5L) rows_int5[[length(rows_int5) + 1L]] <- rec
      }
    }
  }

  calc_sel <- function(rows, level_name) {
    if (length(rows) == 0L) return(NULL)
    dt <- data.table::rbindlist(rows, fill = TRUE)
    if (nrow(dt) == 0L) return(NULL)
    nz <- dt[true_nz == TRUE]; zr <- dt[true_nz == FALSE]
    tpr <- if (nrow(nz) > 0) mean(nz$selected, na.rm = TRUE) else NA_real_
    fpr <- if (nrow(zr) > 0) mean(zr$selected, na.rm = TRUE) else NA_real_
    av  <- tryCatch({
      if (!all(is.na(dt$score)) &&
          length(unique(dt$true_nz)) == 2 &&
          requireNamespace("pROC", quietly = TRUE))
        as.numeric(pROC::auc(pROC::roc(as.numeric(dt$true_nz), dt$score,
                                         quiet = TRUE, direction = "<")))
      else NA_real_
    }, error = function(e) NA_real_)
    data.table::data.table(
      level          = level_name,
      n_candidates   = nrow(dt),
      n_true_nz      = nrow(nz),
      n_selected     = sum(dt$selected, na.rm = TRUE),
      tpr            = tpr,
      fpr            = fpr,
      selection_auc  = av)
  }

  pieces <- list(
    calc_sel(rows_main, "main_effect"),
    calc_sel(rows_int2, "interaction_2way"),
    calc_sel(rows_int3, "interaction_3way"),
    calc_sel(rows_int4, "interaction_4way"),
    calc_sel(rows_int5, "interaction_5way"))
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(data.table::data.table())
  data.table::rbindlist(pieces, fill = TRUE)
}


#' Recovery metrics (bias, RMSE, coverage of credible intervals)
#'
#' @param fit A fitted \code{ctbn_fit} (must have SEs -- typically a MAP
#'   fit with \code{compute_se = TRUE}).
#' @param network The true network spec.
#'
#' @return A data.table with rows per parameter and columns
#'   \code{condition}, \code{parameter}, \code{param_type}, \code{true_val},
#'   \code{est_mean}, \code{ci_lo}, \code{ci_hi}, \code{bias},
#'   \code{sq_err}, \code{covered}.
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' DTw <- prepare_wide(DT)
#'
#' # Recovery metrics need standard errors (compute_se = TRUE).
#' fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 1,
#'                  compute_se = TRUE)
#' rec <- compute_recovery_metrics(fit, network = net)
#' head(rec)
#' }
#'
#' @export
compute_recovery_metrics <- function(fit, network) {

  if (!inherits(fit, c("ctbn_map_fit", "ctbn_stan_fit"))) {
    method <- fit$call_args$method %||% class(fit)[1]
    message(sprintf(
      "compute_recovery_metrics(): not implemented for method '%s' ",
      method),
      "(needs per-coefficient SEs / posteriors from 'map' or 'stan').")
    return(data.table::data.table())
  }

  ca      <- fit$call_args
  conds   <- network$conds
  is_map  <- inherits(fit, "ctbn_map_fit")
  results <- list()

  mk <- function(pt, pn, tv, e, v, lo, hi)
    data.table::data.table(
      condition  = NA_character_,
      parameter  = pn,
      param_type = pt,
      true_val   = tv,
      est_mean   = e,
      est_var    = v,
      ci_lo      = lo,
      ci_hi      = hi,
      bias       = e - tv,
      sq_err     = (e - tv)^2,
      covered    = as.numeric(tv >= lo & tv <= hi))

  for (m in conds) {
    if (is_map) {
      mfit <- fit$map_fits[[m]]; if (is.null(mfit)) next
      bh <- mfit$beta_hat; seh <- mfit$se_hat
      gh <- mfit$gamma_hat
      gse <- mfit$se_gamma %||% rep(NA_real_, mfit$n_gamma)
      x_cols   <- mfit$x_cols
      x_orders <- mfit$x_orders
    } else {
      sfit <- fit$stan_fits[[m]]; if (is.null(sfit)) next
      post <- as.data.frame(sfit)
      beta_cols <- grep("^beta\\[", names(post), value = TRUE)
      gamma_cols <- grep("^gamma\\[", names(post), value = TRUE)
      bh  <- colMeans(post[, beta_cols, drop = FALSE])
      seh <- apply(post[, beta_cols, drop = FALSE], 2, stats::sd)
      gh  <- colMeans(post[, gamma_cols, drop = FALSE])
      gse <- apply(post[, gamma_cols, drop = FALSE], 2, stats::sd)
      influencers <- setdiff(conds, m)
      tmp <- as.data.frame(matrix(0, 1, length(influencers)))
      colnames(tmp) <- influencers
      ir <- build_interaction_cols(tmp, influencers, ca$max_order)
      x_cols <- c(influencers, ir$cols)
      x_orders <- c(rep(0L, length(influencers)), ir$orders)
    }

    infl <- setdiff(conds, m)

    # Main effects
    for (idx in seq_along(infl)) {
      j <- infl[idx]; key <- paste(m, j, sep = "|")
      tv <- if (key %in% names(network$beta_main)) network$beta_main[[key]] else 0
      s <- pmax(seh[idx], 0, na.rm = TRUE)
      r <- mk("main_effect", paste0("b_", j), tv,
              bh[idx], s^2, bh[idx] - 1.96 * s, bh[idx] + 1.96 * s)
      r$condition <- m
      results[[length(results) + 1]] <- r
    }

    # Higher-order interactions
    add_int <- function(ord_int, true_list, slot_label) {
      pos_set <- which(x_orders == ord_int - 1L)
      if (length(pos_set) == 0) return(invisible())
      for (pos in pos_set) {
        cn <- x_cols[pos]
        pts <- strsplit(cn, "_x_")[[1]]
        keys <- vapply(.permutations(pts),
                        function(p) paste(c(m, p), collapse = "|"),
                        character(1))
        mk_k <- keys[keys %in% names(true_list)]
        tv <- if (length(mk_k) > 0) true_list[[mk_k[1]]] else 0
        s <- pmax(seh[pos], 0, na.rm = TRUE)
        r <- mk(slot_label, paste0("b", ord_int, "_", cn), tv,
                bh[pos], s^2, bh[pos] - 1.96 * s, bh[pos] + 1.96 * s)
        r$condition <- m
        results[[length(results) + 1]] <<- r
      }
    }

    if (ca$max_order >= 2L) add_int(2L, network$beta_int2, "interaction_2way")
    if (ca$max_order >= 3L) add_int(3L, network$beta_int3, "interaction_3way")
    if (ca$max_order >= 4L) add_int(4L, network$beta_int4, "interaction_4way")
    if (ca$max_order >= 5L) add_int(5L, network$beta_int5, "interaction_5way")
  }

  data.table::rbindlist(results, fill = TRUE)
}


# ── Internal: ensure each condition has a _event column ────────────────────
#
# Wide-format DT_test (the user-facing layout produced by
# prepare_wide()) does not carry _event columns; those are
# computed inside the fitters. Public metric helpers should be able
# to take that same wide DT directly, so we compute _event on the fly.
# The semantics match the fitter: _event = 1 at a row iff the patient
# is currently in state 0 for that condition AND transitions to state 1
# in the next observed row for the same eid.
.ensure_event_cols <- function(DT, conds) {
  if (!data.table::is.data.table(DT)) DT <- data.table::as.data.table(DT)
  needs <- vapply(conds,
                  function(cn) !paste0(cn, "_event") %in% names(DT),
                  logical(1))
  if (!any(needs)) return(DT)
  DT <- data.table::copy(DT)
  data.table::setorder(DT, eid, time_to_event)
  for (cn in conds[needs]) {
    ec <- paste0(cn, "_event")
    DT[, (ec) := as.numeric(
      get(cn) == 0 & data.table::shift(get(cn), type = "lead") == 1
    ), by = eid]
    DT[is.na(get(ec)), (ec) := 0L]
  }
  DT
}


# ── Internal: enumerate permutations ────────────────────────────────────────
.permutations <- function(x) {
  n <- length(x)
  if (n <= 1L) return(list(x))
  out <- list()
  for (i in seq_len(n)) {
    rest <- x[-i]
    for (p in .permutations(rest))
      out[[length(out) + 1]] <- c(x[i], p)
  }
  out
}
