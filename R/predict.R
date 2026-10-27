.forecast_target <- function(formula, m, prev, process, component, engine, vertex.attr, control) {
  tn <- list(F = list(m), P = list(m), prev = list(prev))
  spec <- .target_spec(tn, 1L, process, component, engine, vertex.attr)
  .make_target(formula, spec, control)
}

.forecast_draw <- function(target, eta, process, engine, n) {
  s <- .aux_stats(list(target), eta, networks = TRUE)[[1]]
  .sim_matrix(s$el, target, process, engine, n)
}

#' Forecast networks
#'
#' Simulates future networks \eqn{y_{T+1}, \dots, y_{T+h}} from the posterior
#' predictive distribution of a fitted separable model. For binary networks
#' (a [bstergm()] fit) each step simulates the formation and persistence
#' networks and assembles \eqn{y_t = y^P_t \cup (y^F_t \setminus y_{t-1})}.
#' For signed networks, a [bstergm()] fit of the interactions and a
#' [bstergm_sign()] fit of the signs are both needed: the interactions are
#' simulated first, then the signs of their formation and persistence
#' networks, and the next network is assembled as
#' \deqn{x_t = x^P_t \cup (x^F_t \setminus x_{t-1}), \qquad
#'   z_t = \{z^P_t \cup (z^F_t \setminus z_{t-1})\} \cap x_t,}{x_t = x^P_t | (x^F_t & !x_{t-1}), z_t = (z^P_t | (z^F_t & !z_{t-1})) & x_t,}
#' where \eqn{z} are the positive ties; all other interactions are negative.
#'
#' Unlike the auxiliary simulations of the sampler, which start at observed
#' networks, these simulations start at the previous network (new
#' interactions starting as negative), so they need to be longer; see
#' `visits`.
#'
#' @param object a [bstergm()] or [bstergm_sign()] fit with both components
#'   (`"F"` and `"P"`).
#' @param interaction for signed networks, the other fit: the two fits may be
#'   given in either order. `NULL` (default) forecasts binary networks from a
#'   [bstergm()] fit.
#' @param horizon number of future time points.
#' @param nsim number of simulated trajectories.
#' @param start the network to start from (a signed or binary adjacency
#'   matrix, or anything accepted by [signed_adjacency()]); by default the
#'   last observed network.
#' @param visits length of each network simulation, in proposals per free
#'   dyad.
#' @param posterior if `TRUE` (default), each trajectory uses a different draw
#'   from the posterior distributions; if `FALSE`, all use the posterior means.
#' @param seed random seed.
#' @param cores number of processes. Each simulation has its own random number
#'   stream, so results do not depend on the number of cores.
#' @param verbose print progress messages.
#' @param ... unused.
#'
#' @details If the fits use a list of formulas (time-varying covariates), the
#'   formula of the last transition is used for all forecasts. Lagged
#'   covariates (`edgecov("lag.active")`, ...) refer to the previous simulated
#'   network.
#'
#' @return An object of class `bstergm_forecast`, a list with
#'   \item{sims}{a list of `nsim` trajectories, each a list of `horizon`
#'     adjacency matrices (signed or binary);}
#'   \item{start}{the starting network;}
#'   \item{signed, horizon, nsim, times}{the type of network, the settings and
#'     the labels of the forecast time points.}
#'   `summary()` gives, for each forecast time point, the mean and 95% interval
#'   of the numbers of ties (by sign) and of triangles (by sign composition).
#' @seealso [simulate.bstergm()] for posterior predictive replicates of the
#'   observed transitions.
#' @examples
#' data(senators70_am)
#' senate <- signed_series(senators70_am[1:3],
#'                         vertex.attr = attr(senators70_am$S99, "senators"))
#' # Very short runs for illustration only
#' fx <- bstergm(~ edges + nodematch("Party"), senate,
#'               iters = 30, burn.in = 15, aux.iters = 500, seed = 1, verbose = FALSE)
#' summary(predict(fx, horizon = 1, nsim = 2, visits = 3, seed = 2, verbose = FALSE))
#' fz <- bstergm_sign(~ edges + gwesf_pos(0.6) + gwese_neg(0.6), senate,
#'                    iters = 30, burn.in = 15, aux.iters = 500, seed = 1, verbose = FALSE)
#' fc <- predict(fz, fx, horizon = 1, nsim = 2, visits = 3, seed = 2, verbose = FALSE)
#' summary(fc)
#' @export
predict.bstergm <- function(object, interaction = NULL, horizon = 1, nsim = 100, start = NULL,
                            visits = 10, posterior = TRUE, seed = NULL,
                            cores = getOption("Bstergm.cores", 1L), verbose = TRUE, ...) {
  fits <- Filter(Negate(is.null), list(object, interaction))
  if (!all(vapply(fits, inherits, TRUE, "bstergm")))
    stop("'object' and 'interaction' must be fits of bstergm() or bstergm_sign().", call. = FALSE)
  procs <- vapply(fits, `[[`, "", "process")
  if (length(fits) == 1) {
    if (procs != "interaction")
      stop("Forecasting signed networks also needs a bstergm() fit of the interactions, ",
           "passed as 'interaction'.", call. = FALSE)
    sfit <- NULL
    xfit <- fits[[1]]
  } else {
    if (!setequal(procs, c("sign", "interaction")))
      stop("predict() needs one bstergm_sign() fit and one bstergm() fit.", call. = FALSE)
    sfit <- fits[[match("sign", procs)]]
    xfit <- fits[[match("interaction", procs)]]
    if (sfit$data$n != xfit$data$n) stop("The two fits must refer to the same nodes.", call. = FALSE)
  }
  for (f in fits) {
    if (!all(c("F", "P") %in% f$components))
      stop("The fits must include the formation and persistence components.", call. = FALSE)
  }
  signed <- !is.null(sfit)
  horizon <- as.integer(horizon)
  nsim <- as.integer(nsim)
  if (horizon < 1 || nsim < 1) stop("'horizon' and 'nsim' must be positive.", call. = FALSE)

  data <- if (signed) sfit$data else xfit$data
  n <- data$n
  dn <- list(data$node.names, data$node.names)
  start <- if (is.null(start)) data$y[[length(data$y)]] else signed_adjacency(start)
  if (!all(dim(start) == n)) stop("'start' must be a network on the ", n, " nodes of the data.", call. = FALSE)
  if (!signed) start <- (start != 0) * 1L
  dimnames(start) <- dn
  vattr <- data$vertex.attr

  last_formula <- function(fit) {
    fl <- .formula_list(fit$formula, fit$n_transitions)
    fl[[length(fl)]]
  }
  fcontrol <- function(fit) {
    ctrl <- fit$control
    ctrl$aux.iters <- NULL
    ctrl$aux.visits <- visits
    ctrl
  }
  pars <- function(fit, U) {
    if (posterior) .pick_draws(fit, U, nsim)
    else matrix(coef(fit)[[U]], nsim, sum(fit$free[[U]]), byrow = TRUE)
  }
  UU <- c(F = "F", P = "P")
  xform <- last_formula(xfit)
  ctrl.x <- fcontrol(xfit)
  draws <- list(interaction = lapply(UU, function(U) pars(xfit, U)))
  if (signed) {
    sform <- last_formula(sfit)
    ctrl.s <- fcontrol(sfit)
    draws$sign <- lapply(UU, function(U) pars(sfit, U))
  }

  if (verbose) message("Forecasting ", horizon, " step", if (horizon > 1) "s", " ahead (", nsim,
                       " trajectories) ...")
  xt0 <- lapply(UU, function(U)
    .forecast_target(xform, start, start, "interaction", U, "standard", vattr, ctrl.x))
  etamaps <- list(interaction = lapply(xt0, function(tg) as.ergm_model(tg$state)$etamap))
  if (signed) {
    st0 <- lapply(UU, function(U)
      .forecast_target(sform, start, start, "sign", U, sfit$engine, vattr, ctrl.s))
    etamaps$sign <- lapply(st0, function(tg) as.ergm_model(tg$state)$etamap)
  }
  eta_of <- function(process, U, k) {
    fit <- if (process == "sign") sfit else xfit
    ergm.eta(.full_theta(fit, U, draws[[process]][[U]][k, ]), etamaps[[process]][[U]])
  }

  one_step <- function(prev, k, xtargets) {
    if (is.null(xtargets)) {
      xtargets <- lapply(UU, function(U)
        .forecast_target(xform, prev, prev, "interaction", U, "standard", vattr, ctrl.x))
    }
    xF <- .forecast_draw(xtargets$F, eta_of("interaction", "F", k), "interaction", "standard", n) > 0
    xP <- .forecast_draw(xtargets$P, eta_of("interaction", "P", k), "interaction", "standard", n) > 0
    x <- xP | (xF & !(prev != 0))
    if (!signed) {
      m <- x * 1L
      dimnames(m) <- dn
      return(m)
    }
    pos <- prev > 0
    zF <- .forecast_draw(.forecast_target(sform, .signed_from_layers(xF & pos, xF, dn), prev, "sign", "F",
                                          sfit$engine, vattr, ctrl.s),
                         eta_of("sign", "F", k), "sign", sfit$engine, n) > 0
    zP <- .forecast_draw(.forecast_target(sform, .signed_from_layers(xP & pos, xP, dn), prev, "sign", "P",
                                          sfit$engine, vattr, ctrl.s),
                         eta_of("sign", "P", k), "sign", sfit$engine, n) > 0
    z <- (zP | (zF & !pos)) & x
    .signed_from_layers(z, x, dn)
  }

  cores <- max(1L, as.integer(cores))
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  jobs <- .sim_jobs(list("all"), nsim, cores, seed)
  res <- .run_jobs(jobs, function(job) {
    lapply(seq_along(job$rows), function(i) {
      k <- job$rows[i]
      .use_stream(job$rngs[[i]])
      traj <- vector("list", horizon)
      prev <- start
      for (h in seq_len(horizon)) {
        traj[[h]] <- one_step(prev, k, if (h == 1L) xt0 else NULL)
        prev <- traj[[h]]
      }
      traj
    })
  }, cores)
  sims <- unlist(res, recursive = FALSE)
  times <- paste0("T+", seq_len(horizon))
  sims <- lapply(sims, function(tr) stats::setNames(tr, times))
  structure(list(sims = sims, start = start, signed = signed, horizon = horizon, nsim = nsim,
                 times = times, posterior = posterior, visits = visits),
            class = "bstergm_forecast")
}

#' @rdname predict.bstergm
#' @param probs quantiles reported by `summary()`.
#' @export
summary.bstergm_forecast <- function(object, probs = c(0.025, 0.975), ...) {
  stat <- if (isFALSE(object$signed)) {
    function(m) c(ties = sum(m != 0) / 2, triangles = sum(diag((m != 0) %*% (m != 0) %*% (m != 0))) / 6)
  } else {
    function(m) c(positive = sum(m > 0) / 2, negative = sum(m < 0) / 2, ties = sum(m != 0) / 2,
                  signed_triad_census(m))
  }
  obs <- stat(object$start)
  out <- lapply(seq_len(object$horizon), function(h) {
    s <- vapply(object$sims, function(tr) stat(tr[[h]]), obs)
    q <- t(apply(s, 1, stats::quantile, probs = probs, names = FALSE))
    colnames(q) <- paste0(100 * probs, "%")
    cbind(start = obs, mean = rowMeans(s), q)
  })
  names(out) <- object$times
  structure(out, class = "summary.bstergm_forecast")
}

#' @export
print.summary.bstergm_forecast <- function(x, digits = 1, ...) {
  for (h in names(x)) {
    cat("Forecast ", h, ":\n", sep = "")
    print(round(x[[h]], digits))
    cat("\n")
  }
  invisible(x)
}

#' @export
print.bstergm_forecast <- function(x, ...) {
  cat("Forecast of ", x$nsim, if (isFALSE(x$signed)) " binary" else " signed", " network trajectories, ",
      x$horizon, " step", if (x$horizon > 1) "s", " ahead (",
      if (x$posterior) "posterior draws" else "posterior means", ")\n\n", sep = "")
  print(summary(x))
  invisible(x)
}
