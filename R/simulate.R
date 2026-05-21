# =============================================================================
# simulate.R  --  Generic CTBN simulator for any number of nodes
# =============================================================================
#
# Generalises the original 10-node multimorbidity simulation (HTN, IHD, MN, ...)
# to support arbitrary network sizes and arbitrary interaction orders up to 5.
# The user provides:
#   * a network specification (list of parents per node, and effect sizes), or
#   * lets the package generate one randomly via make_random_network()
#
# Three high-level entry points:
#   * simulate_ctbn_data()        -- generates a panel-data dataset
#   * make_default_network()      -- the original 10-node multimorbidity network
#   * make_random_network()       -- random DAG with controlled connectivity
# =============================================================================


#' Generate a random network specification for simulation
#'
#' Creates a network of `n_nodes` conditions, where each condition has
#' between `min_parents` and `max_parents` parent influencers chosen at
#' random. All effect sizes (baseline rates, main effects, 2-way to
#' 5-way interactions, covariate gammas) are sampled from sensible
#' defaults so the resulting simulated data mimic real EHR multimorbidity
#' patterns.
#'
#' @param n_nodes Integer, number of conditions in the network (>= 2).
#' @param min_parents Integer, minimum number of parents per condition.
#' @param max_parents Integer, maximum number of parents per condition.
#'   Must be <= n_nodes - 1.
#' @param node_names Character vector of length n_nodes giving custom
#'   node names. Defaults to \code{C1, C2, ...}.
#' @param baseline_range Numeric length-2 vector; baseline rates are
#'   drawn uniformly from this range (then log-transformed).
#' @param main_effect_sd Numeric, standard deviation of main-effect betas.
#' @param int2_prob,int3_prob,int4_prob,int5_prob Numeric in [0,1];
#'   marginal probability that any 2-way / 3-way / 4-way / 5-way
#'   interaction term is non-zero.
#' @param int_effect_sd Numeric, SD of non-zero interaction betas.
#' @param strong_int_prob Probability that a non-zero interaction is a
#'   "strong" effect (mean ~1.0) rather than a "weak" one (mean 0).
#' @param covariate_names Character vector of covariate names (passed
#'   through to gamma generation).
#' @param seed Integer random seed.
#'
#' @return A list with components:
#' \describe{
#'   \item{conds}{Character vector of condition names.}
#'   \item{beta0}{Named numeric vector of log baseline rates.}
#'   \item{graph}{Named list: parents per condition.}
#'   \item{beta_main}{Named numeric vector of main-effect coefficients.}
#'   \item{beta_int2}{Named list of 2-way interaction coefficients.}
#'   \item{beta_int3}{Named list of 3-way coefficients (if any).}
#'   \item{beta_int4}{Named list of 4-way coefficients (if any).}
#'   \item{beta_int5}{Named list of 5-way coefficients (if any).}
#'   \item{gamma}{Named list of per-node covariate coefficients.}
#'   \item{init_prev}{Named numeric vector of initial prevalences.}
#' }
#'
#' @examples
#' net <- make_random_network(n_nodes = 8, max_parents = 3, seed = 1)
#' str(net, max.level = 1)
#'
#' @export
make_random_network <- function(n_nodes        = 10L,
                                 min_parents    = 1L,
                                 max_parents    = 4L,
                                 node_names     = NULL,
                                 baseline_range = c(0.005, 0.20),
                                 main_effect_sd = 0.35,
                                 int2_prob      = 0.40,
                                 int3_prob      = 0.30,
                                 int4_prob      = 0.20,
                                 int5_prob      = 0.10,
                                 int_effect_sd  = 0.30,
                                 strong_int_prob = 0.30,
                                 covariate_names = c("age", "sex_male",
                                                      "smk_current", "smk_former"),
                                 seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  if (n_nodes < 2L)   stop("n_nodes must be at least 2.")
  if (max_parents > n_nodes - 1L)
    stop("max_parents must be <= n_nodes - 1.")

  conds <- if (is.null(node_names)) paste0("C", seq_len(n_nodes)) else node_names
  if (length(conds) != n_nodes) stop("node_names length != n_nodes.")

  # Baseline rates (log scale)
  rates <- stats::runif(n_nodes, baseline_range[1], baseline_range[2])
  beta0 <- setNames(log(rates), conds)

  # Random graph -- choose parents for each node
  graph <- setNames(lapply(conds, function(m) {
    others <- setdiff(conds, m)
    n_pa   <- sample(min_parents:max_parents, 1)
    sample(others, min(n_pa, length(others)))
  }), conds)

  # Main effects
  beta_main <- list()
  for (m in conds) {
    for (j in graph[[m]]) {
      key  <- paste(m, j, sep = "|")
      beta_main[[key]] <- stats::rnorm(1, 0.30, main_effect_sd)
    }
  }

  # Helper: sample an interaction effect (zero w/ given prob, else strong/weak)
  sample_int_effect <- function(p_nonzero, sd_eff = int_effect_sd) {
    if (stats::runif(1) >= p_nonzero) return(0)
    if (stats::runif(1) < strong_int_prob)
      stats::rnorm(1, 1.0, 0.10)
    else
      stats::rnorm(1, 0, sd_eff)
  }

  beta_int2 <- list(); beta_int3 <- list()
  beta_int4 <- list(); beta_int5 <- list()

  for (m in conds) {
    pa <- graph[[m]]
    if (length(pa) >= 2)
      for (pair in utils::combn(pa, 2, simplify = FALSE))
        beta_int2[[paste(c(m, pair), collapse = "|")]] <-
          sample_int_effect(int2_prob)
    if (length(pa) >= 3)
      for (trip in utils::combn(pa, 3, simplify = FALSE))
        beta_int3[[paste(c(m, trip), collapse = "|")]] <-
          sample_int_effect(int3_prob)
    if (length(pa) >= 4)
      for (q in utils::combn(pa, 4, simplify = FALSE))
        beta_int4[[paste(c(m, q), collapse = "|")]] <-
          sample_int_effect(int4_prob)
    if (length(pa) >= 5)
      for (q in utils::combn(pa, 5, simplify = FALSE))
        beta_int5[[paste(c(m, q), collapse = "|")]] <-
          sample_int_effect(int5_prob)
  }

  # Covariate effects (gamma)
  gamma <- setNames(lapply(conds, function(m) {
    setNames(
      stats::rnorm(length(covariate_names), 0.25, 0.10),
      covariate_names)
  }), conds)

  # Initial prevalence (modest)
  init_prev <- setNames(
    pmin(pmax(stats::runif(n_nodes, 0.02, 0.20), 0.01), 0.30),
    conds)

  list(
    conds      = conds,
    beta0      = beta0,
    graph      = graph,
    beta_main  = unlist(beta_main),
    beta_int2  = beta_int2,
    beta_int3  = beta_int3,
    beta_int4  = beta_int4,
    beta_int5  = beta_int5,
    gamma      = gamma,
    init_prev  = init_prev,
    cov_names  = covariate_names
  )
}


#' Default 10-node multimorbidity network
#'
#' Returns the canonical 10-condition network used in the original
#' CTBN multimorbidity studies (HTN, IHD, MN, OA, DM, HL, DVL, RS,
#' DE, CLRD), with hand-tuned baseline rates, parent sets, and effect
#' sizes drawn from UK Biobank empirical analyses.
#'
#' @param seed Integer, random seed for the stochastic 2-way and 3-way
#'   interaction sampling step (to match the original simulation).
#'
#' @return A network specification list compatible with
#'   \code{\link{simulate_ctbn_data}}.
#'
#' @examples
#' net <- make_default_network()
#' net$conds
#'
#' @export
make_default_network <- function(seed = 2026L) {
  conds <- c("HTN","IHD","MN","OA","DM","HL","DVL","RS","DE","CLRD")
  beta0 <- setNames(log(c(
    HTN = 0.150, IHD = 0.080, MN = 0.005, OA = 0.060, DM = 0.050,
    HL  = 0.040, DVL = 0.040, RS = 0.080, DE = 0.070, CLRD = 0.050)), conds)

  graph <- list(
    HTN = c("DM","IHD","RS"), IHD = c("HTN","DM"), MN = c("HTN"),
    OA  = c("DM","HTN"),       DM = c("HTN","IHD"),
    HL  = c("OA","HTN"),       DVL = c("OA","DM"),
    RS  = c("DE","CLRD"),      DE = c("RS","OA"),
    CLRD = c("RS","DM","HTN"))

  beta_main <- c(
    "HTN|DM"   = 0.75, "HTN|IHD" = 0.35, "HTN|RS"   = 0.30,
    "IHD|HTN"  = 0.35, "IHD|DM"  = 0.55,
    "MN|HTN"   = 0.28,
    "OA|DM"    = 0.38, "OA|HTN"  = 0.32,
    "DM|HTN"   = 0.30, "DM|IHD"  = 0.35,
    "HL|OA"    = 0.32, "HL|HTN"  = 0.28,
    "DVL|OA"   = 0.25, "DVL|DM"  = 0.30,
    "RS|DE"    = 0.40, "RS|CLRD" = 0.45,
    "DE|RS"    = 0.42, "DE|OA"   = 0.20,
    "CLRD|RS"  = 0.50, "CLRD|DM" = 0.35, "CLRD|HTN" = 0.28)

  set.seed(seed)
  beta_int2 <- list()
  for (m in names(graph)) {
    pa <- graph[[m]]
    if (length(pa) < 2) next
    for (pair in utils::combn(pa, 2, simplify = FALSE)) {
      key <- paste(m, pair[1], pair[2], sep = "|")
      beta_int2[[key]] <- if (stats::runif(1) < 0.40)
        (if (stats::runif(1) < 0.70) stats::rnorm(1, 0,   0.3)
         else                        stats::rnorm(1, 1.2, 0.1))
      else 0
    }
  }
  set.seed(seed + 99L)
  beta_int3 <- list()
  for (m in names(graph)) {
    pa <- graph[[m]]
    if (length(pa) < 3) next
    for (trip in utils::combn(pa, 3, simplify = FALSE)) {
      key <- paste(c(m, trip), collapse = "|")
      beta_int3[[key]] <- if (stats::runif(1) < 0.30)
        (if (stats::runif(1) < 0.70) stats::rnorm(1, 0,   0.2)
         else                        stats::rnorm(1, 0.9, 0.1))
      else 0
    }
  }

  gamma <- list(
    HTN  = c(age = 0.35, sex_male = 0.20, smk_current = 0.25, smk_former = 0.20),
    IHD  = c(age = 0.30, sex_male = 0.25, smk_current = 0.40, smk_former = 0.25),
    MN   = c(age = 0.25, sex_male = 0.25, smk_current = 0.35, smk_former = 0.20),
    OA   = c(age = 0.45, sex_male = 0.20, smk_current = 0.25, smk_former = 0.25),
    DM   = c(age = 0.25, sex_male = 0.25, smk_current = 0.30, smk_former = 0.22),
    HL   = c(age = 0.40, sex_male = 0.30, smk_current = 0.20, smk_former = 0.20),
    DVL  = c(age = 0.20, sex_male = 0.20, smk_current = 0.20, smk_former = 0.28),
    RS   = c(age = 0.20, sex_male = 0.25, smk_current = 0.25, smk_former = 0.25),
    DE   = c(age = 0.25, sex_male = -0.20, smk_current = 0.20, smk_former = 0.25),
    CLRD = c(age = 0.30, sex_male = 0.20, smk_current = 0.60, smk_former = 0.20))

  init_prev <- c(HTN = 0.16, IHD = 0.08, MN = 0.06, OA = 0.06, DM = 0.05,
                 HL = 0.04, DVL = 0.04, RS = 0.08, DE = 0.07, CLRD = 0.05)

  list(
    conds     = conds,
    beta0     = beta0,
    graph     = graph,
    beta_main = beta_main,
    beta_int2 = beta_int2,
    beta_int3 = beta_int3,
    beta_int4 = list(),
    beta_int5 = list(),
    gamma     = gamma,
    init_prev = init_prev,
    cov_names = c("age", "sex_male", "smk_current", "smk_former")
  )
}


#' Generate true interaction effects for an arbitrary order
#'
#' Helper used internally by \code{\link{make_random_network}} but
#' exposed for users who want to extend an existing network with
#' higher-order interactions.
#'
#' @param graph Named list of parent vectors (as in a network spec).
#' @param order Integer interaction order in [2, 5].
#' @param prob Probability of being non-zero.
#' @param sd Standard deviation of weak non-zero effects.
#' @param strong_prob Probability a non-zero effect is "strong".
#' @param seed Optional integer seed.
#'
#' @return A named list of interaction coefficients keyed
#'   \code{"<target>|<parent1>|<parent2>|..."}.
#'
#' @examples
#' net <- make_random_network(n_nodes = 6, max_parents = 3, seed = 1)
#' int3 <- generate_true_interactions(net$graph, order = 3,
#'                                     prob = 0.5, seed = 1)
#' str(int3, max.level = 1)
#'
#' @export
generate_true_interactions <- function(graph, order = 2L, prob = 0.4,
                                        sd = 0.3, strong_prob = 0.3,
                                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (order < 2L || order > 5L)
    stop("order must be in [2, 5].")
  out <- list()
  for (m in names(graph)) {
    pa <- graph[[m]]
    if (length(pa) < order) next
    for (grp in utils::combn(pa, order, simplify = FALSE)) {
      key <- paste(c(m, grp), collapse = "|")
      out[[key]] <- if (stats::runif(1) < prob)
        (if (stats::runif(1) < strong_prob) stats::rnorm(1, 1.0, 0.10)
         else                                stats::rnorm(1, 0,   sd))
      else 0
    }
  }
  out
}


# ── Internal: compute true rate q for one (target, state, covariates) ───────
.compute_q_true <- function(net, m, X_state, covs, max_order = 3L) {
  lp <- net$beta0[m]
  pa <- net$graph[[m]]
  conds <- net$conds

  # Main effects
  for (j in pa) {
    key <- paste(m, j, sep = "|")
    if (key %in% names(net$beta_main))
      lp <- lp + net$beta_main[[key]] * X_state[j]
  }

  apply_int <- function(int_list, ord) {
    if (length(pa) < ord || max_order < ord) return(0)
    if (length(int_list) == 0) return(0)
    s <- 0
    for (grp in utils::combn(pa, ord, simplify = FALSE)) {
      key <- paste(c(m, grp), collapse = "|")
      if (key %in% names(int_list) && int_list[[key]] != 0) {
        prod_x <- prod(X_state[grp])
        s <- s + int_list[[key]] * prod_x
      }
    }
    s
  }

  lp <- lp + apply_int(net$beta_int2, 2L)
  lp <- lp + apply_int(net$beta_int3, 3L)
  lp <- lp + apply_int(net$beta_int4, 4L)
  lp <- lp + apply_int(net$beta_int5, 5L)

  # Covariate effects
  gamma_m <- net$gamma[[m]]
  if (!is.null(gamma_m)) {
    cv_used <- intersect(names(gamma_m), names(covs))
    if (length(cv_used))
      lp <- lp + sum(gamma_m[cv_used] * covs[cv_used])
  }
  exp(lp)
}


# ── Internal: smoking transition step (default time-varying covariate) ──────
.sim_smoking_step <- function(smk, dt, q01 = 0.04, q12 = 0.06) {
  if (smk == 2L) return(2L)
  if (smk == 0L) return(if (stats::runif(1) < 1 - exp(-q01 * dt)) 1L else 0L)
  if (stats::runif(1) < 1 - exp(-q12 * dt)) 2L else 1L
}

.smk_to_covs <- function(smk)
  c(smk_current = as.numeric(smk == 1L),
    smk_former  = as.numeric(smk == 2L))


# ── Internal: simulate one patient's trajectory ─────────────────────────────
.simulate_patient <- function(eid, net, t_horizon = 10, dt_step = 5,
                                interaction_order = 3L,
                                use_smoking_dynamics = TRUE) {
  conds <- net$conds
  M     <- length(conds)
  X_state <- setNames(stats::rbinom(M, 1, net$init_prev), conds)

  # Default covariate trajectory: age + sex + smoking (mirrors UK Biobank)
  age0     <- log(stats::rnorm(1, 57, 9.39))
  sex_male <- stats::rbinom(1, 1, 0.56)
  smk      <- sample(0:2, 1, prob = c(0.47, 0.11, 0.42))

  records <- list(); t_now <- 0

  while (t_now < t_horizon) {
    t_next   <- min(t_now + dt_step, t_horizon)
    dt       <- t_next - t_now
    smk_covs <- .smk_to_covs(smk)
    covs     <- c(age = age0 + t_now, sex_male = sex_male, smk_covs)
    at_risk  <- conds[X_state == 0]

    if (length(at_risk) == 0) {
      if (use_smoking_dynamics) smk <- .sim_smoking_step(smk, dt)
      t_now <- t_next; next
    }

    q_vec  <- setNames(
      vapply(at_risk, function(m) .compute_q_true(net, m, X_state, covs,
                                                    max_order = interaction_order),
              numeric(1)),
      at_risk)
    Lambda <- sum(q_vec)
    t_cand <- stats::rexp(1, max(Lambda, 1e-10))

    if (t_cand < dt) {
      m_ev <- if (length(at_risk) == 1L) at_risk
              else sample(at_risk, 1, prob = q_vec / Lambda)
      t_ev <- t_now + t_cand
      rec_row <- data.table::data.table(
        eid           = eid,
        time_to_event = t_now,
        dt            = t_ev - t_now)
      for (cn in conds) data.table::set(rec_row, j = cn, value = X_state[cn])
      rec_row$age         <- age0 + t_now
      rec_row$sex_male    <- sex_male
      rec_row$smk_current <- smk_covs["smk_current"]
      rec_row$smk_former  <- smk_covs["smk_former"]
      rec_row$event_cond  <- m_ev
      records[[length(records) + 1]] <- rec_row
      X_state[m_ev] <- 1L
      t_now <- t_ev
      if (all(X_state == 1)) break
    } else {
      rec_row <- data.table::data.table(
        eid           = eid,
        time_to_event = t_now,
        dt            = dt)
      for (cn in conds) data.table::set(rec_row, j = cn, value = X_state[cn])
      rec_row$age         <- age0 + t_now
      rec_row$sex_male    <- sex_male
      rec_row$smk_current <- smk_covs["smk_current"]
      rec_row$smk_former  <- smk_covs["smk_former"]
      rec_row$event_cond  <- NA_character_
      records[[length(records) + 1]] <- rec_row
      if (use_smoking_dynamics) smk <- .sim_smoking_step(smk, dt)
      t_now <- t_next
    }
  }
  if (length(records) == 0) return(NULL)
  data.table::rbindlist(records)
}


#' Simulate a CTBN multimorbidity dataset
#'
#' Generates a panel-data dataset under a user-specified network. Each
#' patient's trajectory of conditions and time-varying covariates is
#' simulated forward in time using exact rate-based event sampling.
#'
#' This generic generator works with any number of nodes (>= 2) and any
#' interaction order from 1 (main effects only) up to 5 (5-way
#' interactions), generalising the original 10-node multimorbidity
#' simulation. The user specifies the desired interaction order via the
#' \code{interaction_order} argument.
#'
#' @param network A network specification as returned by
#'   \code{\link{make_default_network}} or \code{\link{make_random_network}}.
#' @param n_patients Integer, number of patients to simulate.
#' @param t_horizon Numeric, total follow-up duration (years).
#' @param dt_step Numeric, panel grid step (years).
#' @param interaction_order Integer in [1, 5]. Controls which interaction
#'   terms enter the true rate function:
#'   \itemize{
#'     \item 1 -- main effects only
#'     \item 2 -- main + 2-way
#'     \item 3 -- main + 2-way + 3-way (default)
#'     \item 4 -- up to 4-way
#'     \item 5 -- up to 5-way
#'   }
#' @param use_smoking_dynamics Logical, simulate the default smoking
#'   transition process (passes as a time-varying covariate).
#' @param seed Optional integer random seed.
#' @param parallel Logical, parallelise patient simulation across cores.
#' @param n_cores Integer, number of cores when \code{parallel = TRUE}.
#'
#' @return A long-format data.table with one row per at-risk interval
#'   per patient. Columns: \code{eid}, \code{time_to_event}, \code{dt},
#'   one column per condition, plus \code{age}, \code{sex_male},
#'   \code{smk_current}, \code{smk_former}, and \code{event_cond}
#'   (target of the next event, or NA for a censored interval).
#'
#' @examples
#' \dontrun{
#' # Use the default 10-node network with 2-way interactions
#' net <- make_default_network()
#' DT  <- simulate_ctbn_data(net, n_patients = 100,
#'                            interaction_order = 2, seed = 1)
#'
#' # Random 6-node network with 3-way interactions
#' net2 <- make_random_network(n_nodes = 6, max_parents = 3, seed = 1)
#' DT2  <- simulate_ctbn_data(net2, n_patients = 200,
#'                             interaction_order = 3, seed = 1)
#' }
#'
#' @export
simulate_ctbn_data <- function(network,
                                n_patients         = 1000L,
                                t_horizon          = 10,
                                dt_step            = 5,
                                interaction_order  = 3L,
                                use_smoking_dynamics = TRUE,
                                seed               = NULL,
                                parallel           = FALSE,
                                n_cores            = 2L) {

  if (interaction_order < 1L || interaction_order > 5L)
    stop("interaction_order must be in [1, 5].")
  if (!is.null(seed)) set.seed(seed)

  worker <- function(i) {
    .simulate_patient(eid = i, net = network,
                       t_horizon = t_horizon,
                       dt_step = dt_step,
                       interaction_order = interaction_order,
                       use_smoking_dynamics = use_smoking_dynamics)
  }

  if (parallel && n_cores > 1) {
    if (.Platform$OS.type == "unix") {
      out <- parallel::mclapply(seq_len(n_patients), worker,
                                  mc.cores = n_cores,
                                  mc.preschedule = TRUE)
    } else {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl,
        varlist = c(".simulate_patient", ".compute_q_true",
                    ".sim_smoking_step", ".smk_to_covs", "network",
                    "t_horizon", "dt_step", "interaction_order",
                    "use_smoking_dynamics"),
        envir = environment())
      parallel::clusterEvalQ(cl, library(data.table))
      out <- parallel::parLapply(cl, seq_len(n_patients), worker)
    }
  } else {
    out <- lapply(seq_len(n_patients), worker)
  }

  data.table::rbindlist(out, fill = TRUE)
}
