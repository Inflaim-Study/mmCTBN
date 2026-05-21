# =============================================================================
# zzz.R  --  Package globals, internal helpers, NULL operator
# =============================================================================

# Suppress R CMD check NOTEs for data.table NSE columns
utils::globalVariables(c(
  ".", "..CONDS", "..cols", "..x_cols",
  "eid", "time_to_event", "dt", "fold", "event_cond",
  "true_nz", "selected", "score", "param_type",
  "covered_95", "true_nonzero", "bias", "sq_err",
  "selected_g", "marker", "t_event", "t_ev",
  "poisson_ll", "brier", "tdauc", "eval_time",
  "n_test_rows", "n_test_pats", "model", "target",
  "fpr", "tpr", "selection_auc", "level",
  # plot.R bare-name aesthetic mappings
  "logRR", "RR", "name", "influencer", "value",
  "pair_label", "mean_tdauc", "se_tdauc",
  # fit_classic / fit_cph data.table NSE references
  "..conds", ".pa_key", ".ev10", "t_target", "t_parent", "t_last",
  "status", "time",
  # compute_pred_metrics per-patient survival table
  ".row_lp", ".row_y", "T_obs", "marker"
))

#' Null-coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Internal: identify condition columns in a wide-format data.table
#'
#' Distinguishes condition columns from reserved columns and
#' covariates. Used internally by ctbn_fit() and helpers.
#'
#' @keywords internal
#' @noRd
.find_condition_cols <- function(DT_wide,
                                  fixed_covs = character(0),
                                  time_varying_covs = character(0)) {
  all_covs <- c(fixed_covs, time_varying_covs)
  reserved <- c("eid", "time_to_event", "dt", "fold", all_covs)
  ev_cols  <- grep("_event$", names(DT_wide), value = TRUE)
  setdiff(names(DT_wide), c(reserved, ev_cols))
}

#' Internal: convert columns to numeric safely
#' @keywords internal
#' @noRd
.to_num <- function(v) {
  if (is.null(v)) return(NULL)
  if (is.logical(v)) return(as.numeric(v))
  if (is.factor(v) || is.character(v)) return(as.numeric(as.factor(v)) - 1)
  as.numeric(v)
}

#' Format a duration in seconds as h/m/s string
#' @keywords internal
#' @noRd
.format_duration <- function(secs) {
  secs <- as.numeric(secs)
  if (is.na(secs)) return("NA")
  if (secs < 60)   return(sprintf("%.0fs", secs))
  if (secs < 3600) return(sprintf("%.0fm %.0fs", secs %/% 60, secs %% 60))
  sprintf("%.0fh %.0fm", secs %/% 3600, (secs %% 3600) %/% 60)
}

#' Make a simple progress bar string
#' @keywords internal
#' @noRd
.make_progress_bar <- function(current, total, width = 40) {
  frac   <- min(current / max(total, 1), 1)
  filled <- round(frac * width)
  bar    <- paste0(strrep("#", filled), strrep("-", width - filled))
  sprintf("[%s] %3.0f%%  %d/%d", bar, frac * 100, current, total)
}

# Package onLoad: detect optional dependencies once
.mmctbn_env <- new.env(parent = emptyenv())
.mmctbn_env$has_rstan        <- FALSE
.mmctbn_env$has_future       <- FALSE
.mmctbn_env$has_ggplot2      <- FALSE
.mmctbn_env$stan_model_cache <- list()

.onLoad <- function(libname, pkgname) {
  .mmctbn_env$has_rstan   <- requireNamespace("rstan",         quietly = TRUE)
  .mmctbn_env$has_future  <- requireNamespace("future.apply",  quietly = TRUE)
  .mmctbn_env$has_ggplot2 <- requireNamespace("ggplot2",       quietly = TRUE)
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "mmCTBN: Continuous-Time Bayesian Networks for Multimorbidity (v",
    utils::packageVersion("mmCTBN"), ")\n",
    "  Use ctbn_fit() with method = 'map' (fast) or 'stan' (full Bayes).")
}
