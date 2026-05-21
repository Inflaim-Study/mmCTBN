# =============================================================================
# fit_stan.R  --  Bayesian MCMC fit (rstan) for CTBN models
# =============================================================================

#' Fit a CTBN by full Bayesian MCMC (Stan)
#'
#' Fits one Stan model per target node, using Hamiltonian Monte Carlo
#' (NUTS) to sample from the joint posterior over influencer effects,
#' interactions and covariate coefficients.
#'
#' @inheritParams ctbn_map
#' @param chains Integer, number of MCMC chains (default 4).
#' @param iter Integer, total iterations per chain (default 2000).
#' @param warmup Integer, warmup iterations (default 1000).
#' @param seed Integer, random seed.
#' @param parallel Logical, run nodes in parallel via \pkg{future.apply}.
#'
#' @return An object of class \code{c("ctbn_stan_fit", "ctbn_fit")} with
#'   slots \code{beta_matrix}, \code{pip_matrix}, \code{kappa_matrix},
#'   \code{intensity_matrix}, \code{stan_fits}, \code{call_args}, ...
#'
#' @examples
#' \dontrun{
#' net <- make_default_network()
#' DT  <- simulate_ctbn_data(network = net, n_patients = 1000,
#'                            interaction_order = 2, seed = 1)
#' DTw <- prepare_wide(DT)
#'
#' fit <- ctbn_stan(DTw, prior = "horseshoe", max_order = 2,
#'                   fixed_covs = c("intercept","sex_male"),
#'                   time_varying_covs = "age",
#'                   chains = 2, iter = 1000)
#' summary(fit)
#' }
#'
#' @seealso \code{\link{ctbn_map}} for the fast L-BFGS-B / Laplace
#'   approximation; \code{\link{ctbn_fit}} for the unified entry point.
#' @export
ctbn_stan <- function(DT_wide,
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
                       chains            = 4,
                       iter              = 2000,
                       warmup            = 1000,
                       seed              = 42L,
                       parallel          = FALSE,
                       verbose           = TRUE) {

  if (!requireNamespace("rstan", quietly = TRUE))
    stop("ctbn_stan() requires the 'rstan' package.\n",
         "  Install it via:  install.packages('rstan')", call. = FALSE)

  prior <- match.arg(prior)
  stopifnot(data.table::is.data.table(DT_wide))
  stopifnot(all(c("eid","time_to_event","dt") %in% names(DT_wide)))
  if (max_order < 1L || max_order > 5L)
    stop("max_order must be an integer in [1, 5].")

  # Guard against iter <= warmup, which produces zero post-warmup
  # samples and breaks downstream colMeans / setNames calls.
  if (iter <= warmup) {
    new_warmup <- max(1L, as.integer(floor(iter / 2)))
    warning(sprintf(
      "ctbn_stan(): iter=%d <= warmup=%d gives 0 post-warmup samples. ",
      iter, warmup),
      sprintf("Auto-setting warmup=%d. Pass an explicit warmup= to silence this.",
              new_warmup),
      call. = FALSE)
    warmup <- new_warmup
  }
  if (iter - warmup < 10L)
    warning(sprintf(
      "ctbn_stan(): only %d post-warmup samples per chain; ",
      iter - warmup),
      "posterior summaries will be very noisy.", call. = FALSE)

  all_covs   <- c(fixed_covs, time_varying_covs)
  all_conds  <- .find_condition_cols(DT_wide,
                                      fixed_covs = fixed_covs,
                                      time_varying_covs = time_varying_covs)
  n_cond     <- length(all_conds)
  if (n_cond < 2) stop("Need at least 2 condition columns.")

  target_loop <- if (is.null(target_conditions)) all_conds else {
    bad <- setdiff(target_conditions, all_conds)
    if (length(bad)) stop("target_conditions not found: ",
                           paste(bad, collapse = ", "))
    target_conditions
  }

  if (verbose) {
    message(sprintf("ctbn_stan [%s prior, max_order=%d]: %d targets, %d chains",
                    prior, max_order, length(target_loop), chains))
  }

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

  sm <- .get_stan_model(prior)

  # Per-target worker
  fit_one_target <- function(target) {
    .fit_stan_one_target(
      target = target, DT = DT,
      all_conds = all_conds, all_covs = all_covs,
      max_order = max_order, prior = prior, sm = sm,
      theta = theta, a0 = a0, b0 = b0,
      pi0 = pi0, spike_var = spike_var,
      tau0 = tau0, slab_df = slab_df, slab_scale = slab_scale,
      chains = chains, iter = iter, warmup = warmup,
      seed = seed, verbose = verbose)
  }

  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE))
      stop("parallel=TRUE requires 'future.apply'.")
    fit_results <- future.apply::future_lapply(
      target_loop, fit_one_target,
      future.seed = TRUE,
      future.packages = c("rstan", "data.table"))
  } else {
    fit_results <- lapply(target_loop, fit_one_target)
  }
  names(fit_results) <- target_loop

  # ── Assemble matrices ──────────────────────────────────────────────────────
  mk_mat <- function(fill) matrix(fill, n_cond, n_cond,
                                   dimnames = list(all_conds, all_conds))
  beta_matrix      <- mk_mat(0)
  pip_matrix       <- mk_mat(NA_real_)
  kappa_matrix     <- mk_mat(NA_real_)
  intensity_matrix <- mk_mat(0)
  pip_cov_list     <- setNames(vector("list", n_cond), all_conds)
  stan_fits        <- setNames(vector("list", n_cond), all_conds)
  ref_profiles     <- setNames(vector("list", n_cond), all_conds)

  for (target in target_loop) {
    res <- fit_results[[target]]
    if (is.null(res)) next
    influencers <- setdiff(all_conds, target)

    bm <- res$beta_means
    for (idx in seq_along(influencers)) {
      inf <- influencers[idx]
      beta_matrix[inf, target] <- bm[idx]
    }

    if (prior == "spike_slab") {
      for (idx in seq_along(influencers))
        pip_matrix[influencers[idx], target] <- res$pip_b_means[idx]
      pip_cov_list[[target]] <- res$pip_g_named
    } else if (prior %in% c("lasso", "horseshoe")) {
      for (idx in seq_along(influencers))
        kappa_matrix[influencers[idx], target] <- res$kappa_means[idx]
    }

    # Intensities at reference profile
    for (idx in seq_along(influencers)) {
      inf <- influencers[idx]

      skip <- FALSE
      if (variable_select) {
        if (prior == "spike_slab" &&
            !is.na(pip_matrix[inf, target]) &&
            pip_matrix[inf, target] < pip_threshold) skip <- TRUE
        if (prior %in% c("lasso","horseshoe") &&
            !is.na(kappa_matrix[inf, target]) &&
            kappa_matrix[inf, target] > (1 - pip_threshold)) skip <- TRUE
      }
      if (skip) next

      intensity_matrix[inf, target] <- res$intensity_at_inf[[inf]]
    }

    stan_fits[[target]]    <- res$fit
    ref_profiles[[target]] <- res$ref
  }

  structure(
    list(
      method           = "stan_mcmc",
      prior            = prior,
      beta_matrix      = beta_matrix,
      pip_matrix       = pip_matrix,
      kappa_matrix     = kappa_matrix,
      intensity_matrix = intensity_matrix,
      pip_cov_list     = pip_cov_list,
      stan_fits        = stan_fits,
      models           = stan_fits,
      ref_profiles     = ref_profiles,
      pvalue_matrix    = mk_mat(NA_real_),
      call_args        = list(
        method            = "stan",
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
        chains = chains, iter = iter, warmup = warmup
      )
    ),
    class = c("ctbn_stan_fit", "ctbn_fit")
  )
}


# ── Per-target Stan worker (non-exported) ───────────────────────────────────

# Align a numeric matrix to a target column count by right-padding with
# zeros or trimming the extra columns. Used to reconcile posterior
# sample matrices with the expected design dimensions when rstan drops
# or duplicates parameter blocks (e.g. failed chains).
.pad_or_trim_cols <- function(M, target_ncol) {
  if (!is.matrix(M)) M <- as.matrix(M)
  current <- ncol(M)
  if (current == target_ncol) return(M)
  if (current < target_ncol) {
    pad <- matrix(0, nrow = nrow(M), ncol = target_ncol - current)
    return(cbind(M, pad))
  }
  M[, seq_len(target_ncol), drop = FALSE]
}

.fit_stan_one_target <- function(target, DT, all_conds, all_covs, max_order,
                                  prior, sm,
                                  theta, a0, b0, pi0, spike_var,
                                  tau0, slab_df, slab_scale,
                                  chains, iter, warmup, seed, verbose) {
  ec  <- paste0(target, "_event")
  atr <- DT[get(target) == 0 & dt > 0]
  if (nrow(atr) == 0) {
    warning(sprintf("No at-risk rows for '%s' -- skipping.", target))
    return(NULL)
  }
  df <- as.data.frame(atr)
  n  <- nrow(df)
  influencers <- setdiff(all_conds, target)

  # Build X, including interactions
  ir       <- build_interaction_cols(df, influencers, max_order)
  df       <- ir$data
  x_cols   <- c(influencers, ir$cols)
  x_orders <- c(rep(0L, length(influencers)), ir$orders)

  X_mat <- matrix(0, nrow = n, ncol = length(x_cols))
  for (j in seq_along(x_cols)) X_mat[, j] <- .to_num(df[[x_cols[j]]]) %||% 0
  colnames(X_mat) <- x_cols
  Q_stan <- ncol(X_mat)

  # Build Z
  if (length(all_covs) > 0) {
    z_cols <- all_covs
    Z_mat  <- matrix(0, nrow = n, ncol = length(z_cols))
    for (k in seq_along(z_cols)) Z_mat[, k] <- .to_num(df[[z_cols[k]]]) %||% 0
    colnames(Z_mat) <- z_cols
  } else {
    z_cols <- "__intercept__"
    Z_mat  <- matrix(1, nrow = n, ncol = 1)
    colnames(Z_mat) <- z_cols
  }
  P_stan <- ncol(Z_mat)

  # Build reference profile (mean covariates, all conditions = 0)
  ref <- data.frame(dt = 1)
  for (cond in all_conds) ref[[cond]] <- 0
  for (cv in all_covs) {
    col <- atr[[cv]]
    if (is.null(col)) next
    if (is.factor(col) || is.character(col)) {
      tbl       <- sort(table(col), decreasing = TRUE)
      ref[[cv]] <- names(tbl)[1]
      if (is.factor(col)) ref[[cv]] <- factor(ref[[cv]], levels = levels(col))
    } else {
      ref[[cv]] <- mean(col, na.rm = TRUE)
    }
  }

  # Stan data
  stan_data <- list(
    N_rows     = n,
    Q          = Q_stan,
    P          = P_stan,
    N_obs      = as.integer(df[[ec]]),
    X          = X_mat,
    Z          = Z_mat,
    T          = pmax(df[["dt"]], 1e-10),
    beta_order = as.array(x_orders),
    max_order  = as.integer(max_order),
    theta      = theta
  )
  if (prior == "spike_slab") {
    stan_data$pi0 <- pi0; stan_data$a0 <- a0; stan_data$b0 <- b0
    stan_data$spike_var <- spike_var
  } else if (prior %in% c("structured","lasso")) {
    stan_data$a0 <- a0; stan_data$b0 <- b0
  } else if (prior == "horseshoe") {
    stan_data$tau0 <- tau0; stan_data$slab_df <- slab_df
    stan_data$slab_scale <- slab_scale
  }

  fit <- tryCatch(
    rstan::sampling(sm,
                     data    = stan_data,
                     chains  = chains,
                     iter    = iter,
                     warmup  = warmup,
                     seed    = seed,
                     refresh = if (verbose) max(100L, iter %/% 10L) else 0L),
    error = function(e) {
      warning(sprintf("Stan FAILED for '%s': %s", target, e$message))
      NULL
    })
  if (is.null(fit)) return(NULL)

  post <- tryCatch(as.data.frame(fit),
                    error = function(e) {
                      warning(sprintf(
                        "Could not extract posterior for '%s': %s",
                        target, e$message), call. = FALSE)
                      NULL
                    })
  if (is.null(post) || nrow(post) == 0L || ncol(post) == 0L) {
    warning(sprintf("Empty posterior for '%s' -- skipping.", target),
            call. = FALSE)
    return(NULL)
  }

  # ── Helper: extract columns matching a regex and return colMeans
  # vector of *exactly* the requested length, padding with NA when the
  # parameter is absent (which happens, e.g., when iter == warmup, when
  # rstan drops generated quantities, or when a chain failed to start).
  extract_means <- function(prefix, expected_len) {
    cols <- grep(prefix, names(post), value = TRUE)
    if (!length(cols)) return(rep(NA_real_, expected_len))
    samps <- post[, cols, drop = FALSE]
    if (nrow(samps) == 0L) return(rep(NA_real_, expected_len))
    m <- colMeans(samps, na.rm = TRUE)
    if (length(m) == expected_len) return(unname(m))
    # Length mismatch: pad or truncate to expected, with a soft warning.
    if (length(m) < expected_len) {
      out <- rep(NA_real_, expected_len)
      out[seq_along(m)] <- m
      return(out)
    }
    unname(m[seq_len(expected_len)])
  }

  beta_means  <- extract_means("^beta\\[",  Q_stan)
  gamma_means <- extract_means("^gamma\\[", P_stan)

  pip_b_means  <- rep(NA_real_, Q_stan)
  pip_g_named  <- setNames(rep(NA_real_, P_stan), z_cols)
  kappa_means  <- rep(NA_real_, Q_stan)

  if (prior == "spike_slab") {
    pip_b_means <- extract_means("^pip_beta\\[",  Q_stan)
    pip_g_vec   <- extract_means("^pip_gamma\\[", P_stan)
    pip_g_named <- setNames(pip_g_vec, z_cols)
  } else if (prior %in% c("lasso", "horseshoe")) {
    kappa_means <- extract_means("^kappa\\[", Q_stan)
  }

  # Reconstruct posterior sample matrices (now safe -- post is guaranteed
  # to have rows; missing column blocks yield empty matrices instead of
  # crashing). Used for the predictive intensity at the reference profile.
  beta_cols  <- grep("^beta\\[",  names(post), value = TRUE)
  gamma_cols <- grep("^gamma\\[", names(post), value = TRUE)
  beta_samps <- if (length(beta_cols))
    as.matrix(post[, beta_cols, drop = FALSE]) else
    matrix(0, nrow = nrow(post), ncol = Q_stan)
  gamma_samps <- if (length(gamma_cols))
    as.matrix(post[, gamma_cols, drop = FALSE]) else
    matrix(0, nrow = nrow(post), ncol = P_stan)

  # Align columns to Q_stan / P_stan (right-pad with zeros if Stan
  # dropped some, left-pad if it produced extras).
  if (ncol(beta_samps)  != Q_stan)
    beta_samps  <- .pad_or_trim_cols(beta_samps,  Q_stan)
  if (ncol(gamma_samps) != P_stan)
    gamma_samps <- .pad_or_trim_cols(gamma_samps, P_stan)

  # Intensities at reference (one per influencer)
  intensity_at_inf <- list()
  for (idx in seq_along(influencers)) {
    inf <- influencers[idx]
    ref_inf <- ref
    ref_inf[[inf]] <- 1
    df_new <- as.data.frame(ref_inf)
    int_new <- build_interaction_cols(df_new, influencers, max_order)
    df_new  <- int_new$data
    x_new <- vapply(x_cols, function(cv) {
      v <- df_new[[cv]]; if (is.null(v)) 0 else as.numeric(v)[1]
    }, numeric(1))
    z_new <- if ("__intercept__" %in% z_cols) {
      rep(1, P_stan)
    } else {
      vapply(z_cols, function(cv) {
        v <- df_new[[cv]]; if (is.null(v)) 0 else as.numeric(v)[1]
      }, numeric(1))
    }
    intensity_at_inf[[inf]] <- if (nrow(beta_samps) > 0L) {
      eta_samps <- beta_samps %*% x_new + gamma_samps %*% z_new
      mean(exp(eta_samps), na.rm = TRUE)
    } else {
      # Fallback to plug-in (no samples available)
      exp(sum(beta_means * x_new, na.rm = TRUE) +
            sum(gamma_means * z_new, na.rm = TRUE))
    }
  }

  list(
    fit              = fit,
    beta_means       = beta_means[seq_along(influencers)],
    pip_b_means      = pip_b_means,
    pip_g_named      = pip_g_named,
    kappa_means      = kappa_means[seq_along(influencers)],
    intensity_at_inf = intensity_at_inf,
    ref              = ref,
    x_cols           = x_cols,
    x_orders         = x_orders,
    z_cols           = z_cols
  )
}
