# =============================================================================
# Plotting utilities for mmCTBN
# =============================================================================
#
# All plotting functions are guarded with requireNamespace() because
# {ggplot2}, {igraph}, {ggraph}, and {scales} are in Suggests, not Imports.
# =============================================================================


# ---------- internal helpers -------------------------------------------------

.require_ggplot2 <- function(fn = "this plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop(sprintf("Package {ggplot2} is required for %s. ", fn),
         "Install it via install.packages('ggplot2').", call. = FALSE)
}


# Active-edge table for the network and (optionally) the masked heatmap.
# Honours the prior's selection rule (PIP threshold for spike-slab and
# structured; kappa <= 1 - PIP for lasso / horseshoe / classic-like).
.active_pairs <- function(fit, pip_threshold = NULL) {
  if (is.null(pip_threshold))
    pip_threshold <- fit$call_args$pip_threshold %||% 0.5
  bm <- fit$beta_matrix
  pm <- fit$pip_matrix
  km <- fit$kappa_matrix
  pr <- fit$prior

  is_active <- function(inf, tgt) {
    switch(pr,
      spike_slab = !is.na(pm[inf, tgt]) && pm[inf, tgt] >= pip_threshold,
      structured = !is.na(pm[inf, tgt]) && pm[inf, tgt] >= pip_threshold,
      lasso      = !is.na(km[inf, tgt]) && km[inf, tgt] <= (1 - pip_threshold),
      horseshoe  = !is.na(km[inf, tgt]) && km[inf, tgt] <= (1 - pip_threshold),
      classic    = !is.na(pm[inf, tgt]) && pm[inf, tgt] > 0,
      fctbn      = !is.na(pm[inf, tgt]) && pm[inf, tgt] > 0,
      cph        = !is.na(pm[inf, tgt]) && pm[inf, tgt] > 0,
      FALSE)
  }

  im <- fit$intensity_matrix
  rows <- list()
  for (inf in rownames(bm)) for (tgt in colnames(bm)) {
    if (inf == tgt) next
    if (!is_active(inf, tgt)) next
    rows[[length(rows) + 1L]] <- data.frame(
      from        = inf,
      to          = tgt,
      logRR       = bm[inf, tgt],
      RR          = exp(bm[inf, tgt]),
      pip         = if (!is.null(pm)) pm[inf, tgt] else NA_real_,
      kappa       = if (!is.null(km)) km[inf, tgt] else NA_real_,
      intensity   = if (!is.null(im)) im[inf, tgt] else NA_real_,
      stringsAsFactors = FALSE)
  }
  if (length(rows) == 0L)
    return(data.frame(from = character(0), to = character(0),
                      logRR = numeric(0), RR = numeric(0),
                      pip = numeric(0), kappa = numeric(0),
                      intensity = numeric(0)))
  do.call(rbind, rows)
}


# Selection score normalised to [0, 1] for edge width / alpha aesthetics.
.edge_score <- function(active_df, prior) {
  if (prior %in% c("spike_slab", "structured", "classic", "fctbn", "cph"))
    active_df$pip
  else
    1 - active_df$kappa
}


# ---------- network ----------------------------------------------------------

#' Plot a fitted CTBN as a directed network
#'
#' Draws a directed multi-graph of significant influencer -> target edges.
#' Style mirrors the reference figure used in the package paper:
#'
#' * Edges are coloured by **direction**: blue = excitatory (log RR > 0),
#'   red = inhibitory (log RR < 0). Width and alpha encode the
#'   selection score (PIP for spike-slab / structured, 1 - kappa for
#'   lasso / horseshoe).
#' * Edge labels carry the conditional intensity (events per 1000
#'   person-time) at the reference covariate profile, derived from
#'   \code{fit$intensity_matrix}.
#' * Nodes are placed by Fruchterman-Reingold by default (other igraph
#'   layouts via \code{layout=}).
#'
#' @param fit A \code{ctbn_fit} object.
#' @param pip_threshold Selection threshold (default
#'   \code{fit$call_args$pip_threshold} or 0.5).
#' @param layout Character; igraph/ggraph layout name (default
#'   \code{"fr"} for Fruchterman-Reingold; \code{"circle"} is also nice).
#' @param node_size Node diameter (mm in ggraph units; default 12).
#' @param show_intensity_labels Logical; show conditional-intensity
#'   labels on each edge (default TRUE).
#' @param intensity_scale Multiplier applied to \code{intensity_matrix}
#'   for the edge label (default 1000 -> rate per 1000 person-time).
#' @param seed Integer; fixes the random layout seed.
#'
#' @return A \code{ggplot} object (when \pkg{ggraph} is available) or,
#'   invisibly, the active-edge data.frame.
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 2)
#'
#' # Directed network graph (requires ggraph + igraph).
#' plot_network(fit, layout = "circle")
#' }
#'
#' @export
plot_network <- function(fit,
                         pip_threshold         = NULL,
                         layout                = "fr",
                         node_size             = 12,
                         show_intensity_labels = TRUE,
                         intensity_scale       = 1000,
                         seed                  = 123) {

  .require_ggplot2("plot_network")
  edges <- .active_pairs(fit, pip_threshold)
  if (nrow(edges) == 0L) {
    message("No active edges to plot at the chosen threshold.")
    return(invisible(edges))
  }
  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("ggraph", quietly = TRUE)) {
    message("Install {igraph} and {ggraph} for graphical network plots; ",
            "returning the edge table instead.")
    return(invisible(edges))
  }

  edges$direction <- ifelse(edges$logRR > 0, "excitatory", "inhibitory")
  edges$score     <- .edge_score(edges, fit$prior)
  edges$score     <- pmax(pmin(edges$score, 1), 0)
  edges$intensity <- edges$intensity * intensity_scale
  edges$edge_label <- if (show_intensity_labels)
    sprintf("%.1f", edges$intensity) else ""

  dir_colours <- c(excitatory = "#2166ac", inhibitory = "#d73027")
  edges$edge_colour <- dir_colours[edges$direction]

  all_nodes <- rownames(fit$beta_matrix)
  node_df   <- data.frame(name = all_nodes, stringsAsFactors = FALSE)

  g <- igraph::graph_from_data_frame(d = edges, vertices = node_df,
                                       directed = TRUE)

  # ggraph 2.x interpolates edge_colour along arc segments, which
  # forces continuous scale handling and silently drops manual discrete
  # mappings. Pre-assigning the resolved colour via I() bypasses the
  # scale entirely; we attach a separate invisible scale to render a
  # clean direction legend.
  set.seed(seed)
  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_arc(
      ggplot2::aes(edge_colour = I(.data$edge_colour),
                   edge_width  = .data$score,
                   edge_alpha  = .data$score,
                   label       = .data$edge_label),
      arrow       = grid::arrow(length = grid::unit(3, "mm"),
                                 type = "closed"),
      end_cap     = ggraph::circle(node_size * 0.55, "mm"),
      strength    = 0.25,
      angle_calc  = "along",
      label_dodge = grid::unit(3, "mm"),
      label_size  = 2.8,
      label_colour = "grey20",
      show.legend  = FALSE) +
    # Dummy point layer to build a direction colour legend.
    ggplot2::geom_point(
      data    = data.frame(direction = c("excitatory", "inhibitory"),
                            x = NA_real_, y = NA_real_,
                            stringsAsFactors = FALSE),
      ggplot2::aes(x = .data$x, y = .data$y, colour = .data$direction),
      size = 4, na.rm = TRUE, inherit.aes = FALSE) +
    ggplot2::scale_colour_manual(
      values = dir_colours, name = "Direction") +
    ggraph::scale_edge_width(range = c(0.5, 2.8), guide = "none") +
    ggraph::scale_edge_alpha(range = c(0.45, 1),  guide = "none") +
    ggraph::geom_node_point(size = node_size, alpha = 0.88,
                             colour = "#4d4d4d") +
    ggraph::geom_node_text(ggplot2::aes(label = .data$name),
                            colour = "white",
                            size = 3.2,
                            fontface = "bold") +
    ggplot2::labs(
      title    = "Directed multimorbidity network",
      subtitle = sprintf(
        "Edges retained at %.2f selection threshold (%s prior). %s",
        pip_threshold %||% fit$call_args$pip_threshold %||% 0.5,
        fit$prior,
        if (show_intensity_labels)
          paste0("Edge label = conditional intensity per ",
                  intensity_scale, " person-time.")
        else ""),
      caption  = "Blue = excitatory (log RR > 0). Red = inhibitory (log RR < 0).") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      legend.position = "right")
  p
}


# ---------- shared heatmap helper -------------------------------------------

# Long-format conversion that doesn't rely on the .data pronoun.
.matrix_to_long <- function(M, value_name = "value") {
  rn <- rownames(M); cn <- colnames(M)
  if (is.null(rn)) rn <- paste0("R", seq_len(nrow(M)))
  if (is.null(cn)) cn <- paste0("C", seq_len(ncol(M)))
  out <- data.frame(
    influencer = rep(rn, times = length(cn)),
    target     = rep(cn, each  = length(rn)),
    stringsAsFactors = FALSE)
  out[[value_name]] <- as.vector(M)
  out$influencer <- factor(out$influencer, levels = rn)
  out$target     <- factor(out$target,     levels = cn)
  out
}


# ---------- RR heatmap -------------------------------------------------------

#' Heatmap of pairwise rate ratios
#'
#' Diverging heatmap of the \code{beta_matrix} (or its exponential).
#' Cells that fail the prior's selection rule can be greyed out, and
#' the diagonal is always blanked. Active cells are annotated with
#' their RR (or log RR) value.
#'
#' @param fit A \code{ctbn_fit} object.
#' @param scale One of \code{"RR"} (default) or \code{"logRR"}.
#' @param mask_inactive Logical; grey out cells that fail the
#'   selection rule (default TRUE).
#' @param pip_threshold Selection threshold (default
#'   \code{fit$call_args$pip_threshold} or 0.5).
#' @param show_values Logical; annotate active cells with their value
#'   (default TRUE).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 1)
#'
#' plot_heatmap(fit, scale = "RR")
#' }
#'
#' @export
plot_heatmap <- function(fit,
                         scale         = c("RR", "logRR"),
                         mask_inactive = TRUE,
                         pip_threshold = NULL,
                         show_values   = TRUE) {

  .require_ggplot2("plot_heatmap")
  scale <- match.arg(scale)
  bm    <- fit$beta_matrix
  if (is.null(pip_threshold))
    pip_threshold <- fit$call_args$pip_threshold %||% 0.5

  long <- .matrix_to_long(bm, "logRR")
  long$RR <- exp(long$logRR)
  long$same <- long$influencer == long$target

  if (mask_inactive) {
    active <- .active_pairs(fit, pip_threshold)
    if (nrow(active) > 0L) {
      key <- paste(long$influencer, long$target, sep = "::")
      ak  <- paste(active$from, active$to, sep = "::")
      long$logRR[!(key %in% ak)] <- NA_real_
      long$RR   [!(key %in% ak)] <- NA_real_
    } else {
      long$logRR <- NA_real_
      long$RR    <- NA_real_
    }
  }
  long$logRR[long$same] <- NA_real_
  long$RR   [long$same] <- NA_real_

  long$value    <- if (scale == "RR") long$RR else long$logRR
  midpoint      <- if (scale == "RR") 1 else 0
  fill_lab      <- if (scale == "RR") "RR" else "log RR"
  value_format  <- if (scale == "RR") "%.2f" else "%+.2f"

  # Cap fill range at robust quantiles to keep the colour scale stable
  # when a handful of edges have very large rate ratios.
  finite_vals <- long$value[is.finite(long$value)]
  rng <- if (length(finite_vals))
    range(c(midpoint, finite_vals), na.rm = TRUE) else c(midpoint, midpoint)
  if (scale == "RR") rng[1] <- min(rng[1], 1 / max(rng[2], 1))

  p <- ggplot2::ggplot(long,
                       ggplot2::aes(x = .data$target,
                                    y = .data$influencer,
                                    fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_gradient2(
      low = "#d73027", mid = "white", high = "#2166ac",
      midpoint = midpoint, limits = rng,
      na.value = "grey92", name = fill_lab) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_y_discrete(limits = rev(levels(long$influencer))) +
    ggplot2::labs(x = "Target condition", y = "Influencer condition",
                  title = sprintf("Rate-ratio heatmap (%s prior)", fit$prior),
                  subtitle = if (mask_inactive)
                    sprintf("Greyed cells fail the selection rule (threshold %.2f).",
                             pip_threshold) else NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0),
      panel.grid  = ggplot2::element_blank())

  if (show_values) {
    p <- p + ggplot2::geom_text(
      data = long[!is.na(long$value), , drop = FALSE],
      ggplot2::aes(label = sprintf(value_format, .data$value)),
      size = 2.6, colour = "grey15")
  }
  p
}


# ---------- PIP heatmap ------------------------------------------------------

#' Heatmap of posterior inclusion / shrinkage scores
#'
#' For spike-slab / structured priors visualises the PIP matrix; for
#' lasso / horseshoe visualises \code{1 - kappa} (a shrinkage-based
#' inclusion proxy); for the classic / FCTBN / CPH backends visualises
#' the binary \code{pip_matrix} (1 = retained edge).
#'
#' @param fit A \code{ctbn_fit} object.
#' @param threshold Numeric in \eqn{[0, 1]} used to annotate cells
#'   above the threshold (default 0.5).
#' @param show_values Logical; annotate cells whose score >= threshold
#'   (default TRUE).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500, seed = 1)
#' fit <- ctbn_map(prepare_wide(DT), prior = "spike_slab", max_order = 1,
#'                  variable_select = TRUE)
#'
#' plot_pip_heatmap(fit, threshold = 0.5)
#' }
#'
#' @export
plot_pip_heatmap <- function(fit, threshold = 0.5, show_values = TRUE) {

  .require_ggplot2("plot_pip_heatmap")
  pr <- fit$prior
  pm <- fit$pip_matrix
  km <- fit$kappa_matrix

  if (pr %in% c("spike_slab", "structured", "classic", "fctbn", "cph") &&
      !is.null(pm)) {
    M     <- pm
    label <- if (pr == "spike_slab") "PIP" else
              if (pr %in% c("classic", "fctbn", "cph")) "Edge indicator"
              else "pseudo-PIP"
  } else if (pr %in% c("lasso", "horseshoe") && !is.null(km)) {
    M     <- 1 - km
    label <- "1 - kappa"
  } else {
    stop("plot_pip_heatmap(): no inclusion / shrinkage matrix available ",
         "for prior '", pr, "'.")
  }

  long <- .matrix_to_long(M, "value")
  long$value[long$influencer == long$target] <- NA_real_

  p <- ggplot2::ggplot(long,
                       ggplot2::aes(x = .data$target,
                                    y = .data$influencer,
                                    fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_gradient2(
      low = "white", mid = "#9ecae1", high = "#08519c",
      midpoint = 0.3, limits = c(0, 1), na.value = "grey92",
      name = label) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_y_discrete(limits = rev(levels(long$influencer))) +
    ggplot2::labs(x = "Target condition", y = "Influencer condition",
                  title = sprintf("Inclusion heatmap (%s)", label),
                  subtitle = sprintf(
                    "Bold white labels: score >= %.2f.", threshold)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0),
      panel.grid  = ggplot2::element_blank())

  if (show_values) {
    annot <- long[!is.na(long$value) & long$value >= threshold, , drop = FALSE]
    if (nrow(annot)) {
      p <- p + ggplot2::geom_text(
        data = annot,
        ggplot2::aes(label = sprintf("%.2f", .data$value)),
        size = 2.6, colour = "white", fontface = "bold")
    }
  }
  p
}


# ---------- synergy forest ---------------------------------------------------

#' Forest plot of synergistic interaction effects
#'
#' Plots \code{beta_jk} (or the joint log-RR) for the top
#' \code{top_n} interactions by absolute synergistic excess from
#' \code{\link{compute_interaction_effects}}.
#'
#' @param synergy_dt A data.table from \code{compute_interaction_effects}.
#' @param top_n Integer, number of effects to display.
#' @param x One of \code{"beta_jk"} (interaction coef, default) or
#'   \code{"joint_log_RR"} (sum of mains + interaction).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 500,
#'                            interaction_order = 2, seed = 1)
#' DTw <- prepare_wide(DT)
#' fit <- ctbn_map(DTw, prior = "spike_slab", max_order = 2)
#'
#' syn <- compute_interaction_effects(fit, DTw, tau = 5)
#' plot_synergy_forest(syn, top_n = 15)
#' }
#'
#' @export
plot_synergy_forest <- function(synergy_dt, top_n = 25,
                                x = c("beta_jk", "joint_log_RR")) {

  .require_ggplot2("plot_synergy_forest")
  x <- match.arg(x)
  if (!nrow(synergy_dt))
    stop("plot_synergy_forest(): empty synergy_dt.")

  d <- data.table::copy(synergy_dt)
  if (x == "joint_log_RR" && !"joint_log_RR" %in% names(d))
    d$joint_log_RR <- d$beta_j + d$beta_k + d$beta_jk
  d$abs_excess <- abs(d$delta_F %||% d$beta_jk)
  d <- d[order(-d$abs_excess), ]
  if (nrow(d) > top_n) d <- d[seq_len(top_n), ]

  d$pair_label <- paste0(d$condition_j, " x ",
                          d$condition_k, " -> ",
                          d$target_condition)
  d$pair_label <- factor(d$pair_label, levels = rev(d$pair_label))

  ggplot2::ggplot(d, ggplot2::aes(x = .data[[x]], y = .data$pair_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey40") +
    ggplot2::geom_point(size = 2.4, colour = "firebrick") +
    ggplot2::labs(x = if (x == "beta_jk") "Interaction coefficient (log scale)"
                       else "Joint log-RR",
                  y = NULL,
                  title = "Top synergistic interactions") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}


# ---------- tdAUC ------------------------------------------------------------

#' Plot time-dependent AUC from cross-validation
#'
#' Faceted line plot of mean tdAUC vs eval_time, with a 95% across-fold
#' ribbon, for each (model, target).
#'
#' @param cv_results A data.table from \code{\link{ctbn_cv}}.
#' @param targets Optional character subset to display.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' net <- make_random_network(n_nodes = 6, seed = 1)
#' DT  <- simulate_ctbn_data(network = net, n_patients = 400, seed = 1)
#' DTw <- prepare_wide(DT)
#'
#' fit_fns <- list(map = function(DT, ...) ctbn_map(DT, max_order = 1, ...))
#' cv      <- ctbn_cv(DTw, fit_fns = fit_fns, k_folds = 3,
#'                     eval_times = c(1, 3, 5), seed = 1)
#'
#' plot_tdauc(cv)
#' }
#'
#' @export
plot_tdauc <- function(cv_results, targets = NULL) {

  .require_ggplot2("plot_tdauc")
  summ   <- summarise_cv(cv_results)
  auc_dt <- summ$tdauc
  if (!is.null(targets)) auc_dt <- auc_dt[target %in% targets]
  if (!nrow(auc_dt))
    stop("plot_tdauc(): no time-dependent AUC rows in cv_results.")

  ggplot2::ggplot(
    auc_dt,
    ggplot2::aes(x      = .data$eval_time,
                 y      = .data$mean_tdauc,
                 colour = .data$model,
                 fill   = .data$model,
                 group  = .data$model)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$mean_tdauc - 1.96 * .data$se_tdauc,
                   ymax = .data$mean_tdauc + 1.96 * .data$se_tdauc),
      alpha = 0.15, colour = NA) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::facet_wrap(~ target, scales = "free_y") +
    ggplot2::scale_y_continuous(limits = c(0.4, 1)) +
    ggplot2::labs(
      title    = "Time-dependent AUC by model and target (k-fold CV)",
      subtitle = "Ribbon = \u00b11.96 SE across folds",
      x        = "Evaluation time",
      y        = "AUC(t)",
      colour   = "Model",
      fill     = "Model") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}
