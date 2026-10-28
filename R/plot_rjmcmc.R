# Visualisation of RJ-MCMC results for the envidiv model (base R graphics).
#
# Works with the output of rjmcmc_envidiv() and rjmcmc_envidiv_fast(), for a
# single chain or a list of chains (pooled after burn-in).
#
#   s <- rj_water_summary(fit)                 # or list(fit1, fit2, ...)
#   plot_water_posterior(s, phy = phy, ref_stands = literature_stands)
#   plot_water_samples(fit, n = 200)
#   plot_rj_diagnostics(list(fit1, fit2))
#
# Ages are in the units of the sampler (true time); when the time-scale
# factor c was sampled, the tree is drawn with ages divided by the posterior
# median of c.

# ----------------------------------------------------------------- helpers --

#' Wrap a single sampler output in a list
#' @param fits One sampler output or a list of them.
#' @noRd
.as_chain_list <- function(fits) {
  if (!is.null(fits$trace)) list(fits) else fits
}

#' Samples retained after burn-in
#' @param fit Sampler output.
#' @param burnin Fraction discarded.
#' @noRd
.kept_samples <- function(fit, burnin) {
  n <- length(fit$stands)
  keep <- seq.int(floor(burnin * n) + 1, n)
  list(stands = fit$stands[keep], trace = fit$trace[keep, , drop = FALSE])
}

# ------------------------------------------------------------------ summary --

#' Summarise sampled water-level histories
#'
#' Computes, from the output of one or more RJ-MCMC chains (pooled after
#' burn-in), the posterior probability of low water through time, the
#' posterior probability of a drop or a rise in each time bin, and the
#' posterior distribution of the number of lowstands.
#'
#' @param fits Output of [rjmcmc_envidiv_fast()] or [rjmcmc_envidiv()], or a
#'   list of such outputs.
#' @param burnin Fraction of each chain discarded as burn-in.
#' @param grid Ages at which the probability of low water is evaluated; by
#'   default 1001 points from 0 to the crown age.
#' @param bin Width of the time bins (in Myr) for the drop and rise
#'   probabilities.
#' @return A list with elements
#'   \item{low}{`data.frame` with columns `age` and `low_prob`}
#'   \item{drops}{`data.frame` with columns `mid`, `lower`, `upper`,
#'     `drop_prob` and `rise_prob`: the probability that a drop (rise) falls in
#'     each bin}
#'   \item{J}{posterior distribution of the number of lowstands}
#'   \item{T}{crown age in true time}
#'   \item{n_samples}{number of pooled samples}
#'   \item{time_scale}{posterior median of the time-scale factor}
#'   \item{stands, trace}{the pooled samples}
#'   \item{bin}{the bin width}
#' @details Ages are in true time. When the time-scale factor was sampled,
#'   the crown age is divided by its posterior median.
#' @seealso [plot_water_posterior()], [rj_drop_clusters()]
#' @importFrom stats median
#' @importFrom graphics hist
#' @export
rj_water_summary <- function(fits, burnin = 0.2, grid = NULL, bin = 0.01) {
  fits <- .as_chain_list(fits)
  kept <- lapply(fits, .kept_samples, burnin = burnin)
  stands <- unlist(lapply(kept, `[[`, "stands"), recursive = FALSE)
  trace <- do.call(rbind, lapply(kept, `[[`, "trace"))
  T <- max(sapply(fits, `[[`, "T"))
  ts <- if (!is.null(trace$time_scale)) stats::median(trace$time_scale) else 1
  Tt <- T / ts
  S <- length(stands)
  rise <- unlist(lapply(stands, `[[`, "rise"))
  drop <- unlist(lapply(stands, `[[`, "drop"))
  if (is.null(grid)) grid <- seq(0, Tt, length.out = 1001)
  # P(low): lowstands within one sample do not overlap, so the number of
  # samples in low water at age a is #(rise < a) - #(drop < a)
  rs <- sort(rise); ds <- sort(drop)
  n_low <- findInterval(grid, rs, left.open = TRUE) - findInterval(grid, ds)
  low <- data.frame(age = grid, low_prob = n_low / S)
  # drop and rise probabilities per bin
  br <- seq(0, ceiling(Tt / bin) * bin, by = bin)
  dc <- if (length(drop)) graphics::hist(drop, breaks = br, plot = FALSE)$counts else rep(0, length(br) - 1)
  rc <- if (length(rise)) graphics::hist(rise, breaks = br, plot = FALSE)$counts else rep(0, length(br) - 1)
  drops <- data.frame(mid = (br[-1] + br[-length(br)]) / 2, lower = br[-length(br)], upper = br[-1],
                      drop_prob = dc / S, rise_prob = rc / S)
  J <- table(factor(trace$J, levels = 0:max(trace$J))) / nrow(trace)
  list(low = low, drops = drops, J = J, T = Tt, n_samples = S, time_scale = ts,
       stands = stands, trace = trace, bin = bin)
}

#' Credible drop times
#'
#' Merges adjacent time bins with a high posterior drop probability into
#' clusters, each representing one inferred drop of the lake level.
#'
#' @param summary Output of [rj_water_summary()].
#' @param threshold Bins with a drop probability of at least `threshold` are
#'   merged into clusters.
#' @param min_prob Clusters are reported if the posterior probability of at
#'   least one drop in them is at least `min_prob`.
#' @return A `data.frame` with one row per cluster and columns `from`, `to`
#'   (age range), `prob` (posterior probability of at least one drop in the
#'   cluster) and `mean_age` (posterior mean age of those drops).
#' @export
rj_drop_clusters <- function(summary, threshold = 0.05, min_prob = 0.2) {
  d <- summary$drops
  hot <- d$drop_prob >= threshold
  if (!any(hot)) return(data.frame())
  runs <- rle(hot); ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1
  out <- lapply(which(runs$values), function(k) {
    lo <- d$lower[starts[k]]; hi <- d$upper[ends[k]]
    in_cl <- lapply(summary$stands, function(s) s$drop[s$drop >= lo & s$drop < hi])
    p <- mean(lengths(in_cl) > 0)
    data.frame(from = lo, to = hi, prob = p, mean_age = mean(unlist(in_cl)))
  })
  out <- do.call(rbind, out)
  out[out$prob >= min_prob, , drop = FALSE]
}

# --------------------------------------------------------------------- plots --

#' Age axis for plots drawn in time since the crown
#' @param T Crown age.
#' @param side Axis side.
#' @param label Draw labels.
#' @param xlim Age range shown.
#' @noRd
.age_axis <- function(T, side = 1, label = TRUE, xlim = c(T, 0)) {
  at <- pretty(sort(xlim), n = 6)
  at <- at[at >= min(xlim) & at <= max(xlim) & at <= T]
  graphics::axis(side, at = T - at, labels = if (label) format(round(at, 6)) else FALSE)
}

#' Draw reference lowstands as bands
#' @param ref_stands `data.frame` with columns `rise` and `drop`.
#' @param T Crown age.
#' @param ylim Vertical extent of the bands.
#' @param col Fill colour.
#' @noRd
.draw_ref <- function(ref_stands, T, ylim, col) {
  if (is.null(ref_stands) || nrow(ref_stands) == 0) return(invisible())
  graphics::rect(T - ref_stands$drop, ylim[1], T - ref_stands$rise, ylim[2], col = col, border = NA)
}

#' Plot the posterior water-level history
#'
#' Draws, on one age axis, the phylogeny (optional) with the posterior
#' probability of low water as background shading, the posterior probability
#' of low water through time, and the posterior probability of a drop per
#' time bin with the node ages as a rug.
#'
#' @param summary Output of [rj_water_summary()].
#' @param phy The phylogeny (class `phylo`), drawn on top and aligned with the
#'   age axis; `NULL` to omit. When the time-scale factor was sampled, the
#'   tree is rescaled by its posterior median.
#' @param ref_stands `data.frame` with columns `rise` and `drop` of reference
#'   lowstands (for example geological reconstructions or, for simulated
#'   data, the truth), drawn as grey bands with dashed lines at the drops.
#' @param xlim Age range to show, for example `c(1.5, 0)`; by default the
#'   whole tree.
#' @param show_nodes If `TRUE` and `phy` is given, node ages are shown as a
#'   rug below the drop panel.
#' @param col_low,col_drop,col_ref Colours of the low-water probability, the
#'   drop probability and the reference lowstands.
#' @param main Title.
#' @return `summary`, invisibly.
#' @examples
#' \dontrun{
#' s <- rj_water_summary(list(fit1, fit2))
#' plot_water_posterior(s, phy = phy, ref_stands = lowstands)
#' plot_water_posterior(s, xlim = c(1.5, 0))
#' }
#' @importFrom grDevices adjustcolor grey
#' @export
plot_water_posterior <- function(summary, phy = NULL, ref_stands = NULL, xlim = NULL,
                                 show_nodes = TRUE, col_low = "#2c7fb8",
                                 col_drop = "#d95f0e", col_ref = grey(0.85),
                                 main = NULL) {
  T <- summary$T
  if (is.null(xlim)) xlim <- c(T, 0)
  xl <- T - xlim                               # plotting coordinate: time since crown
  n_panels <- if (is.null(phy)) 2 else 3
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  heights <- if (n_panels == 3) c(2.2, 1, 1.2) else c(1, 1.2)
  graphics::layout(matrix(seq_len(n_panels), ncol = 1), heights = heights)
  graphics::par(oma = c(3.5, 0, if (is.null(main)) 0.5 else 2, 0), mgp = c(2.2, 0.6, 0), las = 1)

  node_ages <- NULL
  if (!is.null(phy)) {
    a <- ape::node.depth.edgelength(phy)
    ph <- phy
    ph$edge.length <- ph$edge.length / summary$time_scale     # to true time
    node_ages <- (max(a) - a)[(ape::Ntip(phy) + 1):(ape::Ntip(phy) + phy$Nnode)] / summary$time_scale
    graphics::par(mar = c(0.3, 4.5, 0.5, 1))
    Tphy <- max(a) / summary$time_scale
    # plot.phylo draws the root at x = 0 and the tips at x = Tphy, which is
    # the same coordinate as T - age used in all panels
    yl <- c(0, ape::Ntip(phy) + 1)
    graphics::plot(NA, xlim = xl, ylim = yl, axes = FALSE, xlab = "", ylab = "")
    usr <- graphics::par("usr")
    # low-water posterior as background shading
    lw <- summary$low
    for (i in seq_len(nrow(lw) - 1)) if (lw$low_prob[i] > 0.02)
      graphics::rect(T - lw$age[i + 1], usr[3], T - lw$age[i], usr[4],
                     col = grDevices::adjustcolor(col_low, alpha.f = 0.35 * lw$low_prob[i]), border = NA)
    graphics::par(new = TRUE)
    ape::plot.phylo(ph, show.tip.label = FALSE, x.lim = xl - (T - Tphy), y.lim = yl,
                    edge.width = 0.6, no.margin = FALSE)
    graphics::mtext("tree", side = 2, line = 1, las = 0, cex = 0.8)
  }

  # P(low water)
  graphics::par(mar = c(0.3, 4.5, 0.5, 1))
  graphics::plot(NA, xlim = xl, ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "P(low water)")
  .draw_ref(ref_stands, T, c(0, 1), col_ref)
  lw <- summary$low
  graphics::polygon(c(T - lw$age, rev(T - lw$age)), c(lw$low_prob, rep(0, nrow(lw))),
                    col = grDevices::adjustcolor(col_low, 0.5), border = col_low)
  graphics::axis(2, at = c(0, 0.5, 1)); graphics::box(); .age_axis(T, label = FALSE, xlim = xlim)

  # drop probability per bin
  graphics::par(mar = c(0.5, 4.5, 0.3, 1))
  d <- summary$drops
  ymax <- max(0.05, max(d$drop_prob) * 1.08)
  graphics::plot(NA, xlim = xl, ylim = c(0, ymax), axes = FALSE, xlab = "",
                 ylab = sprintf("P(drop) per %g Myr", summary$bin))
  .draw_ref(ref_stands, T, c(0, ymax), col_ref)
  if (!is.null(ref_stands)) graphics::abline(v = T - ref_stands$drop, col = grey(0.55), lty = 2)
  graphics::rect(T - d$upper, 0, T - d$lower, d$drop_prob, col = col_drop, border = NA)
  if (show_nodes && !is.null(node_ages)) graphics::rug(T - node_ages, ticksize = 0.05, col = grey(0.3))
  graphics::axis(2); graphics::box(); .age_axis(T, xlim = xlim)
  graphics::mtext("Age (Ma)", side = 1, line = 2.2, cex = 0.8)
  if (!is.null(main)) graphics::mtext(main, side = 3, outer = TRUE, line = 0.5, font = 2)
  invisible(summary)
}

#' Plot individual sampled water-level histories
#'
#' Draws a subsample of the posterior water-level histories, one row per
#' sample with low water in colour.
#'
#' @param fits Output of [rjmcmc_envidiv_fast()] or [rjmcmc_envidiv()], or a
#'   list of such outputs.
#' @param n Number of samples drawn (evenly spaced over the pooled samples).
#' @param burnin Fraction of each chain discarded as burn-in.
#' @param ref_stands `data.frame` with columns `rise` and `drop` of reference
#'   lowstands, drawn as grey bands.
#' @param xlim Age range to show; by default the whole tree.
#' @param col_low,col_ref Colours of the sampled and the reference lowstands.
#' @param sort_by_J If `TRUE`, rows are sorted by the number of lowstands.
#' @param main Title.
#' @return The drawn water-level histories, invisibly.
#' @export
plot_water_samples <- function(fits, n = 200, burnin = 0.2, ref_stands = NULL, xlim = NULL,
                               col_low = "#2c7fb8", col_ref = grey(0.85), sort_by_J = TRUE,
                               main = "Sampled water histories") {
  s <- rj_water_summary(fits, burnin = burnin)
  T <- s$T
  if (is.null(xlim)) xlim <- c(T, 0)
  idx <- unique(round(seq(1, s$n_samples, length.out = min(n, s$n_samples))))
  st <- s$stands[idx]
  if (sort_by_J) st <- st[order(sapply(st, function(x) length(x$rise)))]
  graphics::plot(NA, xlim = T - xlim, ylim = c(0, length(st) + 1), axes = FALSE,
                 xlab = "Age (Ma)", ylab = "posterior sample", main = main)
  .draw_ref(ref_stands, T, c(0, length(st) + 1), col_ref)
  for (i in seq_along(st)) if (length(st[[i]]$rise))
    graphics::rect(T - st[[i]]$drop, i - 0.45, T - st[[i]]$rise, i + 0.45, col = col_low, border = NA)
  .age_axis(T, xlim = xlim); graphics::axis(2, at = c(1, length(st))); graphics::box()
  invisible(st)
}

#' MCMC diagnostics for the RJ-MCMC samplers
#'
#' Traces of the log-likelihood, the number of lowstands, the rates, `sigma`,
#' `r` and (if sampled) the time-scale factor, one colour per chain, and the
#' posterior distribution of the number of lowstands per chain.
#'
#' @param fits Output of [rjmcmc_envidiv_fast()] or [rjmcmc_envidiv()], or a
#'   list of such outputs (one per chain).
#' @param burnin Fraction of each chain treated as burn-in (marked in the
#'   traces and excluded from the distribution of the number of lowstands).
#' @return `NULL`, invisibly.
#' @importFrom grDevices hcl.colors
#' @importFrom stats sd
#' @export
plot_rj_diagnostics <- function(fits, burnin = 0.2) {
  fits <- .as_chain_list(fits)
  cols <- grDevices::hcl.colors(max(3, length(fits)), "Dark 3")[seq_along(fits)]
  vars <- c("loglik", "J", "mu", "lambda_H", "lambda_L", "alpha", "sigma", "r")
  has_ts <- any(sapply(fits, function(f) !is.null(f$trace$time_scale) && stats::sd(f$trace$time_scale) > 0))
  if (has_ts) vars <- c(vars, "time_scale")
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  nr <- ceiling((length(vars) + 1) / 2)
  graphics::par(mfrow = c(nr, 2), mar = c(3, 5.5, 1.5, 1), mgp = c(2, 0.6, 0), las = 1)
  for (v in vars) {
    rng <- range(unlist(lapply(fits, function(f) f$trace[[v]])), finite = TRUE)
    log_scale <- v %in% c("mu", "lambda_H", "lambda_L", "alpha", "sigma", "r")
    graphics::plot(NA, xlim = range(unlist(lapply(fits, function(f) f$trace$iter))), ylim = rng,
                   log = if (log_scale) "y" else "", xlab = "iteration", ylab = "", main = "")
    graphics::mtext(v, side = 2, line = 4, las = 0, cex = 0.8)
    for (k in seq_along(fits)) graphics::lines(fits[[k]]$trace$iter, fits[[k]]$trace[[v]], col = cols[k], lwd = 0.6)
    graphics::abline(v = burnin * max(fits[[1]]$trace$iter), lty = 3)
  }
  Js <- lapply(fits, function(f) { k <- .kept_samples(f, burnin); k$trace$J })
  lev <- 0:max(unlist(Js))
  tab <- sapply(Js, function(j) table(factor(j, levels = lev)) / length(j))
  graphics::barplot(t(tab), beside = TRUE, col = cols, names.arg = lev, xlab = "number of lowstands",
                    ylab = "", border = NA)
  graphics::mtext("posterior probability", side = 2, line = 4, las = 0, cex = 0.8)
  if (length(fits) > 1) graphics::legend("topright", fill = cols, legend = paste("chain", seq_along(fits)),
                                         bty = "n", border = NA, cex = 0.8)
  invisible(NULL)
}
