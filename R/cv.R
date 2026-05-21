# =============================================================================
# Cross-validation utilities for mmCTBN
# =============================================================================
#
# Patient-level k-fold CV for any combination of fitted CTBN models. Per-fold
# work is dispatched through `future.apply` so the user's plan() controls
# parallelism (sequential / multisession / multicore). Per-(fold, model,
# target) evaluation is delegated to .eval_cell() so it can also be called
# directly by other parallel backends.
# =============================================================================


#' Construct patient-level fold assignments
#'
#' Splits patients (\code{eid}) — not rows — into \code{k} folds. This
#' avoids leakage of within-patient correlation between train and test.
#'
#' @param DT A data.table with an \code{eid} column.
#' @param k Integer, number of folds (default 5).
#' @param seed Integer, RNG seed (default 42).
#'
#' @return A data.table with columns \code{eid} and \code{fold}.
#'
#' @examples
#' net <- make_random_network(n_nodes = 5, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 60, seed = 1)
#' folds <- make_patient_folds(DT, k = 3, seed = 1)
#' head(folds)
#' table(folds$fold)
#'
#' @export
make_patient_folds <- function(DT, k = 5, seed = 42) {
  if (!"eid" %in% names(DT))
    stop("make_patient_folds(): 'eid' column missing.")
  eids <- unique(DT$eid)
  N    <- length(eids)
  set.seed(seed)
  fold_id <- sample(rep(seq_len(k), length.out = N))
  data.table(eid = eids, fold = fold_id)
}


# ── Per-(fold, model, target) evaluation cell ────────────────────────────────
#
# Returns a data.table of one row per (fold, model, target [, eval_time]).
# Used by .fold_worker() but also exposed (unexported) so users can write
# their own parallel drivers.

.eval_cell <- function(fold_k, model_name, fit, DT_te,
                       all_conds, eval_times) {

  if (is.null(fit)) return(NULL)

  cell_results <- list()

  for (target in all_conds) {

    ec  <- paste0(target, "_event")
    atr <- DT_te[get(target) == 0 & dt > 0]
    if (nrow(atr) == 0) next

    pred <- tryCatch(
      get_lp(fit, atr, target),
      error = function(e) {
        warning(sprintf("get_lp failed [fold=%d, %s, %s]: %s",
                        fold_k, model_name, target, e$message))
        list(lp     = rep(NA_real_, nrow(atr)),
             lambda = rep(NA_real_, nrow(atr)))
      })

    y      <- atr[[ec]]
    dt_    <- atr[["dt"]]
    lp     <- pred$lp
    lambda <- pred$lambda

    mu_hat <- lambda * pmax(dt_, 1e-10)

    poisson_ll <- mean(
      dpois(y, lambda = pmax(mu_hat, 1e-15), log = TRUE),
      na.rm = TRUE)

    p_hat <- 1 - exp(-pmax(mu_hat, 0))
    brier <- mean((y - p_hat)^2, na.rm = TRUE)

    pat_dt <- data.table(
      event   = y,
      t_event = atr$time_to_event,
      marker  = lp + log(pmax(dt_, 1e-10))
    )[is.finite(marker)]

    if (!is.null(eval_times) &&
        nrow(pat_dt) > 1 &&
        length(unique(pat_dt$event)) == 2 &&
        requireNamespace("timeROC", quietly = TRUE)) {

      td <- tryCatch(
        timeROC::timeROC(
          T      = pat_dt$t_event,
          delta  = pat_dt$event,
          marker = pat_dt$marker,
          cause  = 1,
          times  = eval_times,
          iid    = FALSE),
        error = function(e) {
          warning(sprintf("timeROC failed [fold=%d, %s, %s]: %s",
                          fold_k, model_name, target, e$message))
          NULL
        })

      for (t_idx in seq_along(eval_times)) {
        cell_results[[length(cell_results) + 1]] <- data.table(
          fold        = fold_k,
          model       = model_name,
          target      = target,
          poisson_ll  = poisson_ll,
          brier       = brier,
          eval_time   = eval_times[t_idx],
          tdauc       = if (!is.null(td)) td$AUC[t_idx] else NA_real_,
          n_test_rows = nrow(atr),
          n_test_pats = uniqueN(atr$eid))
      }

    } else {
      cell_results[[length(cell_results) + 1]] <- data.table(
        fold        = fold_k,
        model       = model_name,
        target      = target,
        poisson_ll  = poisson_ll,
        brier       = brier,
        eval_time   = NA_real_,
        tdauc       = NA_real_,
        n_test_rows = nrow(atr),
        n_test_pats = uniqueN(atr$eid))
    }
  }

  rbindlist(cell_results)
}


# ── Per-fold worker ──────────────────────────────────────────────────────────

.fold_worker <- function(fold_k, DT_wide, fit_fns, all_conds,
                         eval_times, extra_args, .progress_fn = NULL) {

  DT_tr <- copy(DT_wide[fold != fold_k])[, fold := NULL]
  DT_te <- copy(DT_wide[fold == fold_k])[, fold := NULL]

  setorder(DT_te, eid, time_to_event)
  for (cond in all_conds) {
    ec <- paste0(cond, "_event")
    DT_te[, (ec) := as.numeric(
      get(cond) == 0 & shift(get(cond), type = "lead") == 1
    ), by = eid]
    DT_te[is.na(get(ec)), (ec) := 0L]
  }

  fold_results <- list()

  for (model_name in names(fit_fns)) {

    fit <- tryCatch(
      do.call(fit_fns[[model_name]],
              c(list(DT_tr, target_conditions = all_conds), extra_args)),
      error = function(e) {
        message(sprintf("  [ERROR] Fold %d, model '%s' FAILED: %s",
                        fold_k, model_name, conditionMessage(e)))
        NULL
      })

    cell_dt <- .eval_cell(fold_k, model_name, fit, DT_te,
                          all_conds, eval_times)
    if (!is.null(cell_dt) && nrow(cell_dt) > 0)
      fold_results[[model_name]] <- cell_dt

    if (!is.null(.progress_fn)) .progress_fn()
  }

  rbindlist(fold_results)
}


#' Patient-level k-fold cross-validation for CTBN models
#'
#' Compares any number of CTBN model specifications under a common k-fold
#' patient-level split. Per-fold computation is dispatched via
#' \code{future.apply::future_lapply}, so set the desired plan() before
#' calling, e.g.\ \code{future::plan(multisession, workers = 4)}.
#'
#' @param DT_wide Wide-format data.table (one row per at-risk interval
#'   per patient). Must contain \code{eid}, \code{time_to_event},
#'   \code{dt}, the condition columns and any covariates.
#' @param fit_fns Named list of fitting functions. Each must take a
#'   data.table as its first argument and accept
#'   \code{target_conditions} as a named argument. Most commonly these
#'   are wrappers around \code{ctbn_map} and/or \code{ctbn_stan}.
#' @param fixed_covs Character vector of fixed covariate names.
#' @param time_varying_covs Character vector of time-varying covariate
#'   names.
#' @param k_folds Integer, number of folds (default 3).
#' @param targets Optional character subset of conditions to evaluate.
#' @param seed Integer RNG seed for fold assignment and worker streams.
#' @param eval_times Optional numeric vector of evaluation times for
#'   time-dependent AUC (requires \pkg{timeROC}).
#' @param .progress Logical; if TRUE and \pkg{progressr} is installed,
#'   show progress.
#' @param ... Additional arguments forwarded to every fit function.
#'
#' @return A long-format data.table with columns: \code{fold, model,
#'   target, poisson_ll, brier, eval_time, tdauc, n_test_rows,
#'   n_test_pats}.
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 400,
#'                            interaction_order = 2, seed = 1)
#' DTw <- prepare_wide(DT)
#'
#' # Compare a MAP fit against a classical CTBN under 3-fold CV.
#' # Each fit function takes the data as its first argument and must
#' # accept `target_conditions`; ctbn_cv() injects it per fold.
#' fit_fns <- list(
#'   map     = function(DT, ...) ctbn_map(DT, prior = "spike_slab",
#'                                          max_order = 2, ...),
#'   classic = function(DT, ...) ctbn_classic(DT, max_parents = 3, ...))
#'
#' future::plan(future::sequential)
#' cv <- ctbn_cv(DTw, fit_fns = fit_fns, k_folds = 3,
#'                eval_times = c(2, 5), seed = 1)
#' head(cv)
#' summarise_cv(cv)
#' }
#'
#' @export
ctbn_cv <- function(DT_wide,
                    fit_fns,
                    fixed_covs        = character(0),
                    time_varying_covs = character(0),
                    k_folds           = 3,
                    targets           = NULL,
                    seed              = 42,
                    eval_times        = NULL,
                    .progress         = FALSE,
                    ...) {

  dots <- list(...)
  if ("target_conditions" %in% names(dots))
    stop(paste0(
      "Do not pass 'target_conditions' via '...' or inside fit_fns lambdas.\n",
      "ctbn_cv() injects target_conditions automatically from 'targets'."))

  if (!requireNamespace("future.apply", quietly = TRUE))
    stop("ctbn_cv() requires the {future.apply} package. ",
         "Install it or run sequentially.")

  fold_dt <- make_patient_folds(DT_wide, k = k_folds, seed = seed)
  DT_wide <- merge(copy(DT_wide), fold_dt, by = "eid")

  all_covs       <- c(fixed_covs, time_varying_covs)
  reserved       <- c("eid", "time_to_event", "dt", "fold", all_covs)
  all_conds_full <- setdiff(names(DT_wide), reserved)

  if (is.null(targets)) {
    all_conds <- all_conds_full
  } else {
    missing_targets <- setdiff(targets, all_conds_full)
    if (length(missing_targets))
      stop(sprintf("targets not found in DT_wide: %s",
                   paste(missing_targets, collapse = ", ")))
    all_conds <- targets
  }

  p_fn <- NULL
  if (.progress) {
    if (!requireNamespace("progressr", quietly = TRUE)) {
      warning(".progress = TRUE requires the {progressr} package; ignoring.")
    } else {
      p    <- progressr::progressor(steps = k_folds * length(fit_fns))
      p_fn <- function() p()
    }
  }

  fold_list <- future.apply::future_lapply(
    X              = seq_len(k_folds),
    FUN            = function(fold_k) {
      .fold_worker(
        fold_k       = fold_k,
        DT_wide      = DT_wide,
        fit_fns      = fit_fns,
        all_conds    = all_conds,
        eval_times   = eval_times,
        extra_args   = dots,
        .progress_fn = p_fn)
    },
    future.seed     = seed,
    future.packages = c("data.table", "mmCTBN"))

  rbindlist(fold_list)
}


#' Summarise cross-validation results
#'
#' Aggregates fold-level metrics from \code{\link{ctbn_cv}}.
#'
#' @param cv_results A data.table from \code{ctbn_cv}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{scalar}{data.table with mean/SE of \code{poisson_ll} and
#'     \code{brier} per (model, target).}
#'   \item{tdauc}{data.table with mean/SE of time-dependent AUC per
#'     (model, target, eval_time).}
#' }
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 400, seed = 1)
#' DTw <- prepare_wide(DT)
#'
#' fit_fns <- list(map = function(DT, ...) ctbn_map(DT, max_order = 1, ...))
#' cv      <- ctbn_cv(DTw, fit_fns = fit_fns, k_folds = 3,
#'                     eval_times = c(2, 5), seed = 1)
#'
#' summ <- summarise_cv(cv)
#' summ$scalar   # mean / SE of Poisson log-lik and Brier per model+target
#' summ$tdauc    # mean / SE of time-dependent AUC
#' }
#'
#' @export
summarise_cv <- function(cv_results) {

  scalar_dt <- unique(cv_results[, .(fold, model, target, poisson_ll, brier,
                                     n_test_rows, n_test_pats)])

  scalar_summary <- scalar_dt[, .(
    mean_poisson_ll = mean(poisson_ll, na.rm = TRUE),
    se_poisson_ll   = sd(poisson_ll,   na.rm = TRUE) / sqrt(.N),
    mean_brier      = mean(brier,      na.rm = TRUE),
    se_brier        = sd(brier,        na.rm = TRUE) / sqrt(.N),
    n_folds         = .N
  ), by = .(model, target)][order(target, -mean_poisson_ll)]

  auc_summary <- cv_results[!is.na(eval_time) & !is.na(tdauc), .(
    mean_tdauc = mean(tdauc, na.rm = TRUE),
    se_tdauc   = sd(tdauc,   na.rm = TRUE) / sqrt(.N),
    n_folds    = .N
  ), by = .(model, target, eval_time)][order(target, eval_time, model)]

  list(scalar = scalar_summary, tdauc = auc_summary)
}
