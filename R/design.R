# =============================================================================
# design.R  --  Design-matrix construction for CTBN models
# =============================================================================
#
# The CTBN log-rate for target node m takes the form:
#
#   eta_m = sum_j  beta_jm  * X_j(t)            [main effects, order p=0]
#         + sum_{j<k}  beta_jkm * X_j*X_k       [2-way, p=1]
#         + sum_{j<k<l} beta_jklm * X_j*X_k*X_l [3-way, p=2]
#         + ... up to user-specified max_order
#         + sum_r gamma_rm * Z_r(t)             [covariates]
#
# build_interaction_cols() expands the columns; build_design_matrix() builds
# the full design matrix Phi_m = [X_inf | interactions | Z] for one target.
# =============================================================================


#' Add interaction columns to a data.frame
#'
#' Constructs all pairwise, three-way, ..., k-way interaction columns
#' from a set of base variables and appends them to the data.frame
#' (or data.table). Each new column is named by joining the
#' parent variable names with `"_x_"`.
#'
#' @param data A data.frame or data.table containing the base variables.
#' @param vars Character vector of base variable names to interact.
#' @param k Integer, maximum interaction order (1 = main effects only,
#'   2 = include 2-way, ..., up to length(vars)).
#'
#' @return A list with elements:
#' \describe{
#'   \item{data}{The augmented data.frame including the new columns.}
#'   \item{cols}{Character vector of new interaction column names.}
#'   \item{orders}{Integer vector of interaction orders (1 = 2-way,
#'     2 = 3-way, ...).}
#' }
#'
#' @details The maximum allowed interaction order in `mmCTBN` is 5 (i.e.
#'   `k` may range from 1 to 5). Higher orders are computationally
#'   prohibitive and statistically unidentifiable in typical EHR data.
#'
#' @examples
#' df <- data.frame(A = c(1, 0, 1), B = c(1, 1, 0), C = c(0, 1, 1))
#' res <- build_interaction_cols(df, c("A", "B", "C"), k = 3)
#' res$cols
#'
#' @export
build_interaction_cols <- function(data, vars, k) {
  if (length(vars) < 2 || k < 2)
    return(list(data = data, cols = character(0), orders = integer(0)))
  if (k > 5)
    stop("Interaction order > 5 is not allowed (k must be in 1..5).",
         call. = FALSE)
  k_eff  <- min(k, length(vars))
  cols   <- character(0)
  orders <- integer(0)
  for (ord in 2:k_eff) {
    for (grp in utils::combn(vars, ord, simplify = FALSE)) {
      cn         <- paste(grp, collapse = "_x_")
      data[[cn]] <- Reduce(`*`, lapply(grp, function(v) as.numeric(data[[v]])))
      cols       <- c(cols,   cn)
      orders     <- c(orders, ord - 1L)
    }
  }
  list(data = data, cols = cols, orders = orders)
}


#' Build the per-node design matrix Phi_m
#'
#' Assembles the full design matrix
#' \eqn{\Phi_m = [X_{-m} \mid \text{interactions} \mid Z]}
#' for a single target node \eqn{m}, ready for use in MAP/MCMC fitting.
#'
#' @param X_all Numeric matrix of all condition indicators
#'   (rows = at-risk intervals, columns = condition names).
#' @param Z Numeric matrix of covariates
#'   (rows = at-risk intervals, columns = covariate names).
#' @param cond_names Character vector of all condition names.
#' @param target Character, name of the target condition.
#' @param max_order Integer, maximum interaction order
#'   (1 = main effects only; up to 5 allowed).
#'
#' @return A list with components:
#' \describe{
#'   \item{Phi}{Numeric matrix; full design matrix for target.}
#'   \item{x_cols}{Character vector of column names for the X-block.}
#'   \item{x_orders}{Integer vector of interaction orders for X columns
#'     (0 = main effect, 1 = 2-way, 2 = 3-way, ...).}
#'   \item{n_beta}{Number of X columns (condition + interaction terms).}
#'   \item{n_gamma}{Number of Z columns (covariates).}
#' }
#'
#' @examples
#' # Two conditions plus one covariate, target = "B".
#' X <- matrix(c(1, 0, 1,
#'               0, 1, 1), ncol = 2,
#'             dimnames = list(NULL, c("A", "B")))
#' Z <- matrix(c(0.2, -0.5, 1.1), ncol = 1, dimnames = list(NULL, "age"))
#' dm <- build_design_matrix(X, Z, cond_names = c("A", "B"),
#'                            target = "B", max_order = 1)
#' dm$x_cols
#' dm$Phi
#'
#' @export
build_design_matrix <- function(X_all, Z, cond_names, target, max_order) {
  if (max_order < 1L || max_order > 5L)
    stop("max_order must be an integer in [1, 5].", call. = FALSE)

  influencers <- setdiff(cond_names, target)
  X_inf       <- X_all[, influencers, drop = FALSE]
  x_cols      <- influencers
  x_ord       <- rep(0L, length(influencers))

  if (max_order >= 2L) {
    for (ord in 2L:max_order) {
      grps <- utils::combn(influencers, ord, simplify = FALSE)
      for (grp in grps) {
        cn       <- paste(grp, collapse = "_x_")
        col_vals <- Reduce(`*`, lapply(grp, function(v) X_all[, v]))
        X_inf    <- cbind(X_inf, col_vals)
        x_cols   <- c(x_cols, cn)
        x_ord    <- c(x_ord,  ord - 1L)
      }
    }
  }
  colnames(X_inf) <- x_cols
  list(Phi      = cbind(X_inf, Z),
       x_cols   = x_cols,
       x_orders = x_ord,
       n_beta   = length(x_cols),
       n_gamma  = ncol(Z))
}


#' Prepare a long-format CTBN dataset into wide format
#'
#' Strips the \code{event_cond} bookkeeping column (used during data
#' generation) and ensures all condition columns are integers. The
#' result is suitable for input to \code{\link{ctbn_fit}}.
#'
#' @param DT A data.table in long (interval) format produced by
#'   \code{\link{simulate_ctbn_data}} or compiled from EHR data.
#' @param cond_names Character vector of condition column names. If
#'   \code{NULL}, all integer/0-1 columns excluding the reserved fields
#'   are taken to be conditions.
#'
#' @return A data.table in the canonical wide format expected by
#'   \code{ctbn_fit}.
#'
#' @examples
#' net <- make_random_network(n_nodes = 5, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 50, seed = 1)
#' DTw <- prepare_wide(DT)
#' head(DTw)
#'
#' @export
prepare_wide <- function(DT, cond_names = NULL) {
  if (!data.table::is.data.table(DT)) DT <- data.table::as.data.table(DT)
  DT <- data.table::copy(DT)
  if ("event_cond" %in% names(DT)) DT[, event_cond := NULL]
  if (is.null(cond_names)) {
    reserved   <- c("eid", "time_to_event", "dt")
    cand       <- setdiff(names(DT), reserved)
    cond_names <- cand[vapply(cand, function(c) {
      v <- DT[[c]]
      is.numeric(v) && all(v %in% c(0, 1, NA))
    }, logical(1))]
  }
  if (length(cond_names))
    DT[, (cond_names) := lapply(.SD, as.integer), .SDcols = cond_names]
  DT[]
}
