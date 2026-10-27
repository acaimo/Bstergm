#' Diagnose the length of the auxiliary simulations
#'
#' The approximate exchange algorithm assumes that every auxiliary network,
#' simulated by Metropolis-Hastings starting at the observed network, is a
#' draw from the model. If the simulations are too short, the sampled
#' posterior is closer to a contrastive-divergence approximation than to the
#' posterior of the model. This function runs much longer simulations from the
#' observed networks, at the posterior mean or at posterior draws, and records
#' the trajectories of the model statistics. The trajectories should settle
#' within the number of steps used for fitting (the dotted vertical line in
#' the plot). Large drifts, or chains that move towards and freeze at extreme
#' configurations, indicate that `aux.iters` (or `aux.visits`) should be
#' increased, or that the model is near-degenerate at these parameter values.
#'
#' @param object a [bstergm()] or [bstergm_sign()] fit.
#' @param components components to check.
#' @param transitions indices of the transitions to check.
#' @param visits length of the simulations, in proposals per free dyad.
#' @param nsamples number of recorded states along each simulation.
#' @param ndraws number of parameter values: `1` uses the posterior mean,
#'   larger values use draws spread over the posterior sample.
#' @param seed random seed.
#'
#' @return An object of class `aux_diagnostics_bstergm`: for each component a
#'   list over transitions and parameter values of trajectories of the
#'   statistics, with the observed statistics and the number of visits per
#'   free dyad used for fitting.
#' @seealso [bgof()], [control.bstergm()]
#' @examples
#' data(senators70_am)
#' senate <- signed_series(senators70_am[1:3])
#' fit <- bstergm_sign(~ edges + gwesf_pos(0.6) + gwese_neg(0.6), senate,
#'                iters = 100, burn.in = 50, aux.iters = 1000, seed = 1, verbose = FALSE)
#' dg <- aux_diagnostics(fit, visits = 10, nsamples = 20, seed = 1)
#' dg
#' plot(dg)
#' @export
aux_diagnostics <- function(object, components = object$components, transitions = 1,
                            visits = 50, nsamples = 100, ndraws = 1, seed = NULL) {
  if (!inherits(object, "bstergm")) stop("'object' must be a bstergm fit.", call. = FALSE)
  components <- intersect(components, object$components)
  Tn <- object$n_transitions
  transitions <- intersect(as.integer(transitions), seq_len(Tn))
  if (!length(transitions)) stop("No valid transitions selected.", call. = FALSE)
  formulas <- .formula_list(object$formula, Tn)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  streams <- .rng_streams(seed, length(components))
  out <- list(components = components, transitions = transitions, visits = visits,
              process = object$process, results = list())
  for (k in seq_along(components)) {
    U <- components[k]
    targets <- .build_targets(formulas, object$data, object$process, U, object$engine, object$control)
    m <- as.ergm_model(targets[[1]]$state)
    draws <- if (ndraws <= 1) matrix(coef(object)[[U]], 1) else .pick_draws(object, U, ndraws)
    res <- .run_jobs(list(list(rng = streams[[k]])), function(job) {
      lapply(transitions, function(t) {
        tg <- targets[[t]]
        interval <- max(1L, as.integer(round(visits * tg$nfree / nsamples)))
        lapply(seq_len(nrow(draws)), function(r) {
          eta <- ergm.eta(.full_theta(object, U, draws[r, ]), m$etamap)
          z <- ergm_MCMC_slave(tg$state, eta = eta, control = tg$control, verbose = FALSE,
                               burnin = 0L, samplesize = as.integer(nsamples), interval = interval)
          s <- rbind(tg$obs, z$s)
          list(visits = c(0, seq_len(nsamples) * interval / tg$nfree), stats = s)
        })
      })
    })[[1]]
    names(res) <- paste0("t", transitions)
    out$results[[U]] <- list(traj = res, obs = lapply(targets[transitions], `[[`, "obs"),
                             used = object$aux.iters[[U]][transitions] /
                               vapply(targets[transitions], `[[`, 1, "nfree"),
                             names = param_names(m, canonical = TRUE))
  }
  class(out) <- "aux_diagnostics_bstergm"
  out
}

#' @rdname aux_diagnostics
#' @param x an `aux_diagnostics_bstergm` object.
#' @param ... unused.
#' @export
print.aux_diagnostics_bstergm <- function(x, ...) {
  cat("Auxiliary simulations started at the observed networks (", x$visits,
      " proposals per free dyad)\n", sep = "")
  cat("drift = (mean over the second half - value at the fitting length) / SD over the second half\n",
      "        (SD floored at 10% of the SD of the whole trajectory)\n", sep = "")
  for (U in x$components) {
    r <- x$results[[U]]
    for (j in seq_along(r$traj)) {
      tr <- r$traj[[j]][[1]]
      s <- tr$stats
      half <- tr$visits >= max(tr$visits) / 2
      at <- s[which.min(abs(tr$visits - r$used[j])), ]
      sdv <- pmax(apply(s[half, , drop = FALSE], 2, sd), 0.1 * apply(s, 2, sd))
      tab <- cbind(observed = s[1, ], `at fitting length` = at,
                   `second half mean` = colMeans(s[half, , drop = FALSE]),
                   drift = (colMeans(s[half, , drop = FALSE]) - at) / ifelse(sdv > 0, sdv, NA))
      rownames(tab) <- r$names
      cat("\nComponent ", U, ", transition ", names(r$traj)[j], " (fitting used ",
          format(round(r$used[j], 2)), " proposals per free dyad)\n", sep = "")
      print(round(tab, 2))
    }
  }
  invisible(x)
}

#' @rdname aux_diagnostics
#' @export
plot.aux_diagnostics_bstergm <- function(x, ...) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  for (U in x$components) {
    r <- x$results[[U]]
    for (j in seq_along(r$traj)) {
      p <- length(r$names)
      par(mfrow = grDevices::n2mfrow(p), mar = c(3.2, 3.2, 2, 0.6), mgp = c(1.9, 0.6, 0),
          oma = c(0, 0, 2, 0))
      for (q in seq_len(p)) {
        rng <- range(unlist(lapply(r$traj[[j]], function(tr) tr$stats[, q])))
        plot(NA, xlim = c(0, x$visits), ylim = rng, xlab = "proposals per free dyad",
             ylab = "statistic", main = r$names[q], cex.main = 0.8, ...)
        for (tr in r$traj[[j]]) lines(tr$visits, tr$stats[, q], col = adjustcolor("grey20", 0.7))
        abline(h = r$traj[[j]][[1]]$stats[1, q], col = "black", lty = 2)
        abline(v = r$used[j], col = "grey50", lty = 3, lwd = 1.5)
      }
      mtext(paste0(.component_label(U, x$process), ", transition ", names(r$traj)[j]), outer = TRUE, font = 2)
    }
  }
  invisible(x)
}
