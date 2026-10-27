.component_label <- function(U, process) {
  paste0(c(F = "Formation", P = "Persistence")[[U]], " (", U, ") - ",
         if (process == "sign") "sign process" else "binary network")
}

.draws_matrix <- function(object, U) {
  a <- object$theta[[U]]
  m <- matrix(aperm(a, c(1, 3, 2)), ncol = dim(a)[2])
  colnames(m) <- dimnames(a)[[2]]
  m
}

#' @export
print.bstergm <- function(x, digits = 3, ...) {
  cat("Bayesian separable temporal ERGM for signed networks\n")
  cat("  process:", x$process, "| method:", toupper(x$method), "| engine:", x$engine, "\n")
  f <- if (inherits(x$formula, "formula")) x$formula else x$formula[[1]]
  cat("  formula:", paste(deparse(f, width.cutoff = 500L), collapse = " "), "\n")
  cat("  transitions:", x$n_transitions, "| iterations:", x$iters, "| burn-in:", x$burn.in,
      "| thin:", x$thin, "\n")
  for (U in x$components) {
    d <- .draws_matrix(x, U)
    cat("\n", .component_label(U, x$process), ": ", x$n_chains[[U]], " chain",
        if (x$n_chains[[U]] > 1) "s", ", acceptance ", paste(format(round(x$acceptance[[U]], 3)), collapse = " "),
        "\n", sep = "")
    tab <- cbind(Mean = colMeans(d), SD = apply(d, 2, sd))
    print(round(tab, digits))
  }
  invisible(x)
}

#' Summarise a bstergm fit
#'
#' Posterior summaries of each fitted component: mean, standard deviation,
#' naive and time-series standard errors of the mean, quantiles, effective
#' sample size and, with several chains, the potential scale reduction factor
#' of Gelman and Rubin.
#'
#' @param object a [bstergm()] or [bstergm_sign()] fit.
#' @param probs quantiles to report.
#' @param ... unused.
#' @return An object of class `summary.bstergm`: a list with one table per
#'   component plus the acceptance rates.
#' @export
summary.bstergm <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  tabs <- list()
  for (U in object$components) {
    ml <- as.mcmc.list.bstergm(object, components = U, prefix = FALSE)
    d <- as.matrix(ml)
    ess <- tryCatch(coda::effectiveSize(ml), error = function(e) rep(NA_real_, ncol(d)))
    tsse <- apply(d, 2, sd) / sqrt(pmax(ess, 1))
    q <- t(apply(d, 2, quantile, probs = probs))
    tab <- cbind(Mean = colMeans(d), SD = apply(d, 2, sd), `Naive SE` = apply(d, 2, sd) / sqrt(nrow(d)),
                 `Time-series SE` = tsse, q, ESS = ess)
    if (coda::nchain(ml) > 1) {
      rh <- tryCatch(coda::gelman.diag(ml, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1],
                     error = function(e) rep(NA_real_, ncol(d)))
      tab <- cbind(tab, Rhat = rh)
    }
    rownames(tab) <- colnames(d)
    tabs[[U]] <- tab
  }
  out <- list(tables = tabs, acceptance = object$acceptance, process = object$process,
              method = object$method, components = object$components, n_chains = object$n_chains,
              iters = object$iters, burn.in = object$burn.in, thin = object$thin,
              fixed = object$fixed, time = object$time$total)
  class(out) <- "summary.bstergm"
  out
}

#' @export
print.summary.bstergm <- function(x, digits = 3, ...) {
  cat("Posterior summary (", toupper(x$method), ", ", x$iters, " iterations, ", x$burn.in,
      " burn-in, thin ", x$thin, ")\n", sep = "")
  for (U in x$components) {
    cat("\n", .component_label(U, x$process), "\n", sep = "")
    cat("  chains:", x$n_chains[[U]], "| acceptance:",
        paste(format(round(x$acceptance[[U]], 3)), collapse = " "), "\n")
    tab <- x$tables[[U]]
    print(round(tab, digits))
    if (!is.null(x$fixed[[U]]))
      cat("  Fixed (not estimated):", paste(names(x$fixed[[U]]), "=", format(x$fixed[[U]]), collapse = ", "), "\n")
  }
  cat("\nTotal time:", .format_time(x$time), "\n")
  invisible(x)
}

#' Posterior means of a bstergm fit
#'
#' @param object a [bstergm()] or [bstergm_sign()] fit.
#' @param ... unused.
#' @return A list of named vectors of posterior means, one per component.
#' @export
coef.bstergm <- function(object, ...) {
  lapply(object$theta, function(a) apply(a, 2, mean))
}

#' Convert bstergm draws to coda objects
#'
#' @param x a [bstergm()] or [bstergm_sign()] fit.
#' @param components components to include; with both, the parameter names
#'   are prefixed by `F:` and `P:`.
#' @param prefix whether to prefix the parameter names with the component.
#' @param ... unused.
#' @return An [coda::mcmc.list()] with one chain per sampler chain, or, for
#'   `as.mcmc()`, an [coda::mcmc()] object with the chains stacked.
#' @method as.mcmc.list bstergm
#' @export
as.mcmc.list.bstergm <- function(x, components = x$components, prefix = length(components) > 1, ...) {
  components <- intersect(components, x$components)
  H <- unique(vapply(components, function(U) dim(x$theta[[U]])[3], 1L))
  if (length(H) != 1) stop("The components have different numbers of chains; convert them separately.", call. = FALSE)
  start <- x$burn.in + x$thin
  coda::mcmc.list(lapply(seq_len(H), function(h) {
    m <- do.call(cbind, lapply(components, function(U) {
      d <- x$theta[[U]][, , h, drop = FALSE]
      d <- matrix(d, nrow = dim(d)[1], dimnames = list(NULL, dimnames(d)[[2]]))
      if (prefix) colnames(d) <- paste0(U, ":", colnames(d))
      d
    }))
    coda::mcmc(m, start = start, thin = x$thin)
  }))
}

#' @rdname as.mcmc.list.bstergm
#' @method as.mcmc bstergm
#' @export
as.mcmc.bstergm <- function(x, components = x$components, prefix = length(components) > 1, ...) {
  ml <- as.mcmc.list.bstergm(x, components = components, prefix = prefix)
  coda::mcmc(as.matrix(ml))
}

#' Plot the posterior draws of a bstergm fit
#'
#' For each parameter, shows the posterior density, the trace of every chain
#' and the autocorrelation function.
#'
#' @param x a [bstergm()] or [bstergm_sign()] fit.
#' @param components components to plot.
#' @param which any of `"density"`, `"trace"` and `"acf"`.
#' @param lag.max maximum lag of the autocorrelation plots.
#' @param ask whether to ask before each new page (one page per component).
#' @param ... further graphical parameters passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.bstergm <- function(x, components = x$components, which = c("density", "trace", "acf"),
                         lag.max = 50, ask = length(components) > 1 && interactive(), ...) {
  which <- match.arg(which, several.ok = TRUE)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask), add = TRUE)
  }
  for (U in intersect(components, x$components)) {
    a <- x$theta[[U]]
    p <- dim(a)[2]; H <- dim(a)[3]
    par(mfrow = c(p, length(which)), mar = c(3.2, 3.2, 1.6, 0.6), mgp = c(1.9, 0.6, 0),
        oma = c(0, 0, 2, 0))
    cols <- if (H > 1) grDevices::gray.colors(H, start = 0.05, end = 0.6) else "black"
    for (j in seq_len(p)) {
      nm <- dimnames(a)[[2]][j]
      v <- as.vector(a[, j, ])
      if ("density" %in% which) {
        dd <- stats::density(v)
        plot(dd, main = nm, xlab = "", ylab = "Density", ...)
        abline(v = mean(v), lty = 2, col = "grey40")
      }
      if ("trace" %in% which) {
        matplot(seq_len(dim(a)[1]) * x$thin + x$burn.in, matrix(a[, j, ], ncol = H), type = "l", lty = 1,
                col = adjustcolor(cols, 0.8), main = nm, xlab = "Iteration", ylab = "Value", ...)
      }
      if ("acf" %in% which) {
        ac <- stats::acf(a[, j, 1], lag.max = lag.max, plot = FALSE)
        plot(ac$lag[, 1, 1], ac$acf[, 1, 1], type = "h", main = nm, xlab = "Lag", ylab = "ACF",
             ylim = c(min(0, ac$acf), 1), ...)
        abline(h = 0, col = "grey40")
      }
    }
    mtext(.component_label(U, x$process), outer = TRUE, line = 0.5, font = 2)
  }
  invisible(x)
}

#' @importFrom coda as.mcmc.list as.mcmc
#' @export
coda::as.mcmc.list

#' @export
coda::as.mcmc

#' @importFrom ergm.multi Layer
#' @export
ergm.multi::Layer
