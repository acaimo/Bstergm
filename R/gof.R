.full_theta <- function(object, U, draw) {
  th <- object$theta.fixed[[U]]
  th[object$free[[U]]] <- draw
  th
}

.common_prefix <- function(x) {
  if (length(x) == 1) return(x)
  ch <- strsplit(x, "")
  n <- min(lengths(ch))
  k <- 0
  while (k < n && length(unique(vapply(ch, `[`, "", k + 1))) == 1) k <- k + 1
  out <- sub("[^[:alnum:]_)]+$", "", substr(x[1], 1, k))
  if (nzchar(out)) out else x[1]
}

.pick_draws <- function(object, U, nsim) {
  d <- .draws_matrix(object, U)
  idx <- unique(round(seq(1, nrow(d), length.out = min(nsim, nrow(d)))))
  if (length(idx) < nsim) idx <- rep_len(idx, nsim)
  d[idx, , drop = FALSE]
}

#' Posterior predictive goodness of fit
#'
#' Simulates, for posterior draws of the parameters, the formation and
#' persistence networks of every transition (starting from the observed
#' networks, as in the auxiliary step of the sampler) and compares their
#' statistics with the observed ones. Following Caimo and Gollini, the
#' centred cumulative differences
#' \eqn{\sum_t \{s(z^{U\prime}_t) - s(z^U_t)\}} of the model statistics should
#' be centred around zero for a well-fitting model. Additional statistics not
#' in the model (e.g. degree distributions or the signed triangle census) can
#' be checked with `monitor`.
#'
#' @param object a [bstergm()] or [bstergm_sign()] fit.
#' @param nsim number of posterior predictive simulations.
#' @param monitor optional one-sided formula of additional statistics, e.g.
#'   `~ degree(0:10) + degree_neg(0:10) + triangle_signed()` for the sign
#'   process.
#' @param aux.iters number of steps of each simulation; by default five times
#'   the number used for fitting.
#' @param components components to check.
#' @param cores number of processes. Each simulation has its own random number
#'   stream, so results do not depend on the number of cores.
#' @param seed random seed.
#' @param verbose print progress messages.
#'
#' @return An object of class `bgof_bstergm` with, for each component, an
#'   array `sim` (simulations x statistics x transitions) of simulated
#'   statistics, the matrix `obs` of observed statistics and the matrix
#'   `diff` of centred cumulative differences.
#' @seealso [plot.bgof_bstergm()], [bstergm()]
#' @examples
#' data(senators70_am)
#' senate <- signed_series(senators70_am[1:3],
#'                         vertex.attr = attr(senators70_am$S99, "senators"))
#' fit <- bstergm_sign(~ edges + nodematch("Party") + gwesf_pos(0.6) + gwese_pos(0.6),
#'                data = senate, iters = 100, burn.in = 50, aux.iters = 1000, seed = 1,
#'                verbose = FALSE)
#' gf <- bgof(fit, nsim = 10, monitor = ~ triangle_signed(), seed = 2, verbose = FALSE)
#' gf
#' plot(gf)
#' @export
bgof <- function(object, nsim = 100, monitor = NULL, aux.iters = NULL,
                 components = object$components, cores = getOption("Bstergm.cores", 1L),
                 seed = NULL, verbose = TRUE) {
  if (!inherits(object, "bstergm")) stop("'object' must be a bstergm fit.", call. = FALSE)
  components <- intersect(components, object$components)
  Tn <- object$n_transitions
  formulas <- .formula_list(object$formula, Tn)
  if (!is.null(monitor)) {
    if (!inherits(monitor, "formula")) stop("'monitor' must be a one-sided formula.", call. = FALSE)
    if (length(monitor) == 3) monitor <- monitor[-2]
  }
  ctrl <- object$control

  if (verbose) message("Simulating ", nsim, " posterior predictive replicates ...")
  setup <- list()
  for (U in components) {
    ctrlU <- ctrl
    if (is.null(aux.iters)) ctrlU$aux.iters <- NULL else ctrlU$aux.iters <- aux.iters
    targets <- .build_targets(formulas, object$data, object$process, U, object$engine, ctrlU, monitor)
    if (is.null(aux.iters)) for (t in seq_along(targets))
      targets[[t]]$aux.iters <- 5L * object$aux.iters[[U]][t]
    m.model <- ergm_model(formulas[[1]], targets[[1]]$basis)
    m.ext <- as.ergm_model(targets[[1]]$state)
    nmodel <- nparam(m.model, canonical = FALSE)
    ncan <- nparam(m.model, canonical = TRUE)
    groups <- NULL
    if (!is.null(monitor)) {
      m.mon <- ergm_model(monitor, targets[[1]]$basis)
      mterms <- Filter(function(tm) length(tm$coef.names) > 0, m.mon$terms)
      sizes <- vapply(mterms, function(tm) length(tm$coef.names), 1L)
      pos <- cumsum(c(ncan, sizes))
      groups <- lapply(seq_along(mterms), function(k) (pos[k] + 1):pos[k + 1])
      names(groups) <- vapply(mterms, function(tm) .common_prefix(tm$coef.names), "")
    }
    setup[[U]] <- list(targets = targets, etamap = m.ext$etamap, nmodel = nmodel, ncan = ncan,
                       next.par = nparam(m.ext, canonical = FALSE), groups = groups,
                       names = param_names(m.ext, canonical = TRUE),
                       draws = .pick_draws(object, U, nsim))
  }

  cores <- max(1L, as.integer(cores))
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  jobs <- .sim_jobs(components, nsim, cores, seed)

  run_job <- function(job) {
    st <- setup[[job$component]]
    lapply(seq_along(job$rows), function(i) {
      r <- job$rows[i]
      .use_stream(job$rngs[[i]])
      th <- numeric(st$next.par)
      th[seq_len(st$nmodel)] <- .full_theta(object, job$component, st$draws[r, ])
      eta <- ergm.eta(th, st$etamap)
      vapply(.aux_stats(st$targets, eta), `[[`, numeric(length(st$names)), "stats")
    })
  }
  res <- .run_jobs(jobs, run_job, cores)

  out <- list(components = components, process = object$process, nsim = nsim,
              sim = list(), obs = list(), diff = list(), model = list(), groups = list())
  for (U in components) {
    st <- setup[[U]]
    sims <- unlist(res[vapply(jobs, `[[`, "", "component") == U], recursive = FALSE)
    arr <- array(NA_real_, c(nsim, length(st$names), Tn),
                 dimnames = list(NULL, st$names, paste0("t", seq_len(Tn))))
    for (k in seq_len(nsim)) arr[k, , ] <- sims[[k]]
    obs <- .obs_matrix(st$targets, st$names)
    out$sim[[U]] <- arr
    out$obs[[U]] <- obs
    out$diff[[U]] <- apply(sweep(arr, c(2, 3), t(obs)), c(1, 2), sum)
    out$model[[U]] <- seq_len(st$ncan)
    out$groups[[U]] <- st$groups
  }
  class(out) <- "bgof_bstergm"
  out
}

#' @rdname bgof
#' @param x a `bgof_bstergm` object.
#' @param ... unused.
#' @export
print.bgof_bstergm <- function(x, ...) {
  s <- summary.bgof_bstergm(x)
  cat("Posterior predictive goodness of fit (", x$nsim, " simulations)\n", sep = "")
  for (U in x$components) {
    cat("\n", .component_label(U, x$process), "\n", sep = "")
    print(round(s[[U]], 3))
  }
  invisible(x)
}

#' @rdname bgof
#' @return `summary()` returns, for each component, a table with the observed
#'   statistics summed over transitions, the mean and 95% interval of the
#'   simulated ones and a two-sided posterior predictive p-value.
#' @export
summary.bgof_bstergm <- function(object, ...) {
  out <- list()
  for (U in object$components) {
    tot <- apply(object$sim[[U]], c(1, 2), sum)
    obs <- colSums(object$obs[[U]])
    pv <- vapply(seq_along(obs), function(j)
      min(1, 2 * min(mean(tot[, j] >= obs[j]), mean(tot[, j] <= obs[j]))), 0)
    out[[U]] <- cbind(observed = obs, mean = colMeans(tot),
                      `2.5%` = apply(tot, 2, quantile, 0.025), `97.5%` = apply(tot, 2, quantile, 0.975),
                      `p-value` = pv)
  }
  out
}

#' Plot posterior predictive goodness of fit
#'
#' `type = "model"` shows, for each component, boxplots of the centred
#' cumulative differences \eqn{\sum_t \{s(z'_t) - s(z_t)\}} of the model
#' statistics (standardised by their simulated standard deviation when
#' `standardize = TRUE`), which should be centred on the dashed zero line.
#' `type = "monitor"` shows the simulated monitored statistics summed over
#' transitions, with the observed values as black points.
#'
#' @param x a [bgof()] result.
#' @param type `"model"` and/or `"monitor"`.
#' @param components components to plot.
#' @param standardize standardise the centred differences of the model
#'   statistics, so that they share one axis.
#' @param ... further graphical parameters passed to [graphics::boxplot()].
#' @return `x`, invisibly.
#' @export
plot.bgof_bstergm <- function(x, type = c("model", "monitor"), components = x$components,
                              standardize = TRUE, ...) {
  type <- match.arg(type, several.ok = TRUE)
  components <- intersect(components, x$components)
  if (!any(vapply(x$groups[components], length, 1L) > 0)) type <- setdiff(type, "monitor")
  npan <- length(components) * (("model" %in% type) + if ("monitor" %in% type)
    max(vapply(x$groups[components], length, 1L)) else 0)
  op <- par(mfrow = grDevices::n2mfrow(max(npan, 1)), mar = c(7, 4, 2.5, 1), mgp = c(2.2, 0.6, 0))
  on.exit(par(op))
  for (U in components) {
    if ("model" %in% type) {
      d <- x$diff[[U]][, x$model[[U]], drop = FALSE]
      ylab <- "centred cumulative difference"
      if (standardize) {
        s <- apply(apply(x$sim[[U]][, x$model[[U]], , drop = FALSE], c(1, 2), sum), 2, sd)
        d <- sweep(d, 2, ifelse(s > 0, s, 1), "/")
        ylab <- paste(ylab, "(standardised)")
      }
      boxplot(d, las = 2, cex.axis = 0.7, main = paste0(U, ": model statistics"), ylab = ylab,
              col = "grey90", outline = FALSE, ...)
      abline(h = 0, col = "black", lwd = 2, lty = 2)
    }
    if ("monitor" %in% type) {
      tot <- apply(x$sim[[U]], c(1, 2), sum)
      obs <- colSums(x$obs[[U]])
      for (g in names(x$groups[[U]])) {
        idx <- x$groups[[U]][[g]]
        boxplot(tot[, idx, drop = FALSE], las = 2, cex.axis = 0.7, main = paste0(U, ": ", g),
                ylab = "simulated (sum over transitions)", col = "grey90", outline = FALSE,
                names = colnames(tot)[idx], ...)
        points(seq_along(idx), obs[idx], col = "black", pch = 19)
        if (length(idx) > 1) lines(seq_along(idx), obs[idx], col = "black")
      }
    }
  }
  invisible(x)
}
