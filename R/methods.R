# =============================================================================
# methods.R  --  S3 methods for ctbn_fit objects
# =============================================================================

#' Print a CTBN fit
#'
#' @param x A \code{ctbn_fit} object.
#' @param digits Integer, decimal places.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly. Called for the side effect of printing.
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 1)
#'
#' print(fit)
#' }
#'
#' @export
print.ctbn_fit <- function(x, digits = 4, ...) {
  ca     <- x$call_args
  method <- if (!is.null(ca$method)) ca$method else
    if (inherits(x, "ctbn_map_fit"))     "map"     else
    if (inherits(x, "ctbn_stan_fit"))    "stan"    else
    if (inherits(x, "ctbn_classic_fit")) "classic" else
    if (inherits(x, "ctbn_fctbn_fit"))   "fctbn"   else
    if (inherits(x, "ctbn_cph_fit"))     "cph"     else
    "unknown"

  prior_str <- x$prior %||% "n/a"
  cat(sprintf("\n=== mmCTBN Fit (%s%s) ===\n",
              toupper(method),
              if (prior_str %in% c("spike_slab", "structured",
                                    "lasso", "horseshoe"))
                paste0(", ", prior_str, " prior") else ""))
  cat(sprintf("  Conditions      : %d\n", nrow(x$beta_matrix)))
  if (!is.null(ca$max_order))
    cat(sprintf("  Max order       : %d\n", ca$max_order))
  if (!is.null(ca$theta))
    cat(sprintf("  theta (penalty) : %.2f\n", ca$theta))

  switch(prior_str,
    spike_slab = cat(sprintf("  pi0=%.2f | spike_var=%.4f\n",
                              ca$pi0 %||% NA_real_,
                              ca$spike_var %||% NA_real_)),
    horseshoe  = cat(sprintf("  tau0=%.2f | slab_df=%.1f | slab_scale=%.2f\n",
                              ca$tau0 %||% NA_real_,
                              ca$slab_df %||% NA_real_,
                              ca$slab_scale %||% NA_real_)),
    NULL
  )

  if (!is.null(x$convergence)) {
    nc <- sum(x$convergence == 0L, na.rm = TRUE)
    nt <- sum(!is.na(x$convergence))
    cat(sprintf("  Convergence     : %d/%d nodes\n", nc, nt))
  }

  cat("\n-- Beta matrix (influencer -> target) --\n")
  print(round(x$beta_matrix, digits))

  if (prior_str == "spike_slab" && !is.null(x$pip_matrix)) {
    cat("\n-- PIP matrix --\n")
    print(round(x$pip_matrix, digits))
  } else if (prior_str %in% c("lasso", "horseshoe") &&
             !is.null(x$kappa_matrix)) {
    cat("\n-- Kappa matrix (0=signal, 1=shrunk) --\n")
    print(round(x$kappa_matrix, digits))
  } else if (!is.null(x$pip_matrix) && any(is.finite(x$pip_matrix))) {
    cat("\n-- Edge indicator matrix (1 = retained) --\n")
    print(round(x$pip_matrix, digits))
  }

  invisible(x)
}

#' @rdname print.ctbn_fit
#' @export
print.ctbn_map_fit <- function(x, digits = 4, ...) {
  print.ctbn_fit(x, digits = digits, ...)
  if (!is.null(x$se_matrix)) {
    cat("\n-- Laplace SE matrix --\n")
    print(round(x$se_matrix, digits))
  }
  invisible(x)
}

#' @rdname print.ctbn_fit
#' @export
print.ctbn_stan_fit <- function(x, digits = 4, ...) {
  print.ctbn_fit(x, digits = digits, ...)
  ca <- x$call_args
  cat(sprintf("\n  MCMC: %d chains x %d iter (warmup %d)\n",
              ca$chains, ca$iter, ca$warmup))
  invisible(x)
}

#' @rdname print.ctbn_fit
#' @export
print.ctbn_classic_fit <- function(x, digits = 4, ...) {
  print.ctbn_fit(x, digits = digits, ...)
  ca <- x$call_args
  npar <- vapply(x$rates, nrow, integer(1))
  cat(sprintf("\n  Classical CTBN: %d targets, %d total parent configs\n",
              length(x$rates), sum(npar)))
  invisible(x)
}

#' @rdname print.ctbn_fit
#' @export
print.ctbn_fctbn_fit <- function(x, digits = 4, ...) {
  print.ctbn_fit(x, digits = digits, ...)
  n_sel <- sum(x$pip_matrix == 1, na.rm = TRUE)
  n_tot <- sum(!is.na(x$pip_matrix))
  cat(sprintf("\n  FCTBN: %d/%d edges retained by group lasso\n",
              n_sel, n_tot))
  invisible(x)
}

#' @rdname print.ctbn_fit
#' @export
print.ctbn_cph_fit <- function(x, digits = 4, ...) {
  print.ctbn_fit(x, digits = digits, ...)
  n_edges <- sum(vapply(x$beta_cov, length, integer(1)))
  alpha   <- x$call_args$alpha %||% 0.05
  cat(sprintf("\n  CTBN-CPH: %d Cox-PH models (alpha = %.3f)\n",
              n_edges, alpha))
  invisible(x)
}


#' Summary of active influencer pairs
#'
#' @param object A \code{ctbn_fit} object.
#' @param pip_threshold Selection threshold; default uses
#'   \code{object$call_args$pip_threshold}.
#' @param show_covariates Logical, also list active covariates.
#' @param ... Currently unused.
#'
#' @return A list summarising the active influencer pairs, returned
#'   invisibly; printed for inspection.
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 1)
#'
#' summary(fit)
#' }
#'
#' @export
summary.ctbn_fit <- function(object, pip_threshold = NULL,
                              show_covariates = FALSE, ...) {

  if (is.null(pip_threshold))
    pip_threshold <- object$call_args$pip_threshold %||% 0.5

  bm <- object$beta_matrix
  pm <- object$pip_matrix
  km <- object$kappa_matrix
  im <- object$intensity_matrix
  pr <- object$prior

  is_active <- function(inf, tgt) {
    switch(pr,
      spike_slab = !is.na(pm[inf,tgt]) && pm[inf,tgt] >= pip_threshold,
      structured = if (inherits(object, "ctbn_map_fit"))
                     !is.na(pm[inf,tgt]) && pm[inf,tgt] >= pip_threshold
                   else TRUE,
      lasso      = !is.na(km[inf,tgt]) && km[inf,tgt] <= (1 - pip_threshold),
      horseshoe  = !is.na(km[inf,tgt]) && km[inf,tgt] <= (1 - pip_threshold),
      classic    = !is.na(pm[inf,tgt]) && pm[inf,tgt] > 0,
      fctbn      = !is.na(pm[inf,tgt]) && pm[inf,tgt] > 0,
      cph        = !is.na(pm[inf,tgt]) && pm[inf,tgt] > 0,
      FALSE)
  }

  rows <- list()
  for (inf in rownames(bm)) for (tgt in colnames(bm)) {
    if (inf == tgt || !is_active(inf, tgt)) next
    row <- data.frame(
      influencer = inf, target = tgt,
      log_RR     = round(bm[inf,tgt], 4),
      RR         = round(exp(bm[inf,tgt]), 4),
      intensity  = round(im[inf,tgt], 6),
      stringsAsFactors = FALSE)
    if (inherits(object, "ctbn_map_fit") && !is.null(object$se_matrix)) {
      row$SE <- round(object$se_matrix[inf, tgt], 4)
      row$z  <- round(bm[inf,tgt] /
                       pmax(object$se_matrix[inf,tgt], 1e-10), 3)
    }
    if (pr == "spike_slab")
      row$PIP <- round(pm[inf,tgt], 3)
    else if (pr == "horseshoe")
      row$kappa <- round(km[inf,tgt], 3)
    else if (pr == "lasso")
      row$kappa <- round(km[inf,tgt], 3)
    else if (pr == "structured" && !is.null(pm) && !is.na(pm[inf,tgt]))
      row$pseudo_PIP <- round(pm[inf,tgt], 3)
    rows[[length(rows) + 1]] <- row
  }

  sel_lbl <- switch(pr,
    spike_slab = sprintf("PIP >= %.2f", pip_threshold),
    lasso      = sprintf("kappa <= %.2f", 1 - pip_threshold),
    horseshoe  = sprintf("kappa <= %.2f", 1 - pip_threshold),
    structured = "all pairs",
    classic    = "all graph edges",
    fctbn      = "non-zero group lasso coefficients",
    cph        = "all graph edges (covariate-modulated)",
    "")

  cat(sprintf("\n=== Active Influencer Pairs (%s, %s) ===\n", pr, sel_lbl))

  result_inf <- NULL
  if (length(rows) == 0) {
    cat("  None\n")
  } else {
    result_inf <- do.call(rbind, rows)
    sc <- if (pr == "spike_slab") "PIP"
          else if (pr %in% c("lasso","horseshoe")) "kappa"
          else "log_RR"
    if (sc %in% names(result_inf))
      result_inf <- result_inf[order(result_inf[[sc]],
                                       decreasing = (sc == "PIP")), ]
    rownames(result_inf) <- NULL
    print(result_inf)
  }

  result_cov <- NULL
  if (show_covariates && !is.null(object$pip_cov_list)) {
    cov_col_name <- if (pr == "spike_slab") "PIP" else "pseudo_PIP"
    cov_rows <- list()
    for (tgt in names(object$pip_cov_list)) {
      pips <- object$pip_cov_list[[tgt]]
      if (is.null(pips) || length(pips) == 0) next
      for (cv in names(pips)) {
        pv <- pips[[cv]]
        if (!is.na(pv) && pv >= pip_threshold)
          cov_rows[[length(cov_rows) + 1]] <- data.frame(
            covariate = cv, target = tgt,
            stat      = round(pv, 3),
            stringsAsFactors = FALSE)
      }
    }
    if (length(cov_rows)) {
      result_cov           <- do.call(rbind, cov_rows)
      names(result_cov)[3] <- cov_col_name
      result_cov           <- result_cov[order(-result_cov[[cov_col_name]]), ]
      rownames(result_cov) <- NULL
      cat(sprintf("\n=== Active Covariates (%s >= %.2f) ===\n",
                  cov_col_name, pip_threshold))
      print(result_cov)
    }
  }

  invisible(list(influencers = result_inf, covariates = result_cov))
}

#' @rdname summary.ctbn_fit
#' @export
summary.ctbn_map_fit  <- function(object, ...) summary.ctbn_fit(object, ...)
#' @rdname summary.ctbn_fit
#' @export
summary.ctbn_stan_fit <- function(object, ...) summary.ctbn_fit(object, ...)
#' @rdname summary.ctbn_fit
#' @export
summary.ctbn_classic_fit <- function(object, ...) summary.ctbn_fit(object, ...)
#' @rdname summary.ctbn_fit
#' @export
summary.ctbn_fctbn_fit   <- function(object, ...) summary.ctbn_fit(object, ...)
#' @rdname summary.ctbn_fit
#' @export
summary.ctbn_cph_fit     <- function(object, ...) summary.ctbn_fit(object, ...)


#' Extract coefficients
#'
#' @param object A \code{ctbn_fit} object.
#' @param type Which coefficient table to return:
#'   \code{"beta"} (default; condition-condition log-RR),
#'   \code{"intensity"}, \code{"pip"}, \code{"kappa"}, or \code{"se"}.
#' @param ... Currently unused.
#'
#' @return A numeric matrix of the requested coefficient table (may be
#'   \code{NULL} if the fit does not provide it).
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 1)
#'
#' coef(fit)                 # condition-condition log-RR matrix
#' coef(fit, "intensity")    # reference-profile intensities
#' }
#'
#' @export
coef.ctbn_fit <- function(object, type = c("beta", "intensity", "pip",
                                             "kappa", "se"), ...) {
  type <- match.arg(type)
  switch(type,
    beta      = object$beta_matrix,
    intensity = object$intensity_matrix,
    pip       = object$pip_matrix,
    kappa     = object$kappa_matrix,
    se        = object$se_matrix
  )
}


#' Plot a fitted CTBN
#'
#' Default \code{plot} dispatches to \code{\link{plot_network}}.
#'
#' @param x A \code{ctbn_fit} object.
#' @param type One of \code{"network"} (DAG), \code{"heatmap"} (RR
#'   heatmap), \code{"pip"} (PIP heatmap).
#' @param ... Passed to the underlying plot function.
#'
#' @return A \code{ggplot} object from the underlying plot function.
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 1)
#'
#' plot(fit, type = "network")
#' plot(fit, type = "heatmap")
#' }
#'
#' @export
plot.ctbn_fit <- function(x, type = c("network","heatmap","pip"), ...) {
  type <- match.arg(type)
  switch(type,
    network = plot_network(x, ...),
    heatmap = plot_heatmap(x, ...),
    pip     = plot_pip_heatmap(x, ...))
}
