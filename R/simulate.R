.sim_matrix <- function(el, target, process, engine, n) {
  el <- as.matrix(el)[, 1:2, drop = FALSE]
  if (process == "sign") {
    m <- matrix(0L, n, n)
    if (engine == "native") {
      act <- target$basis %n% ".active"
      if (length(act)) m[act] <- m[act[, 2:1, drop = FALSE]] <- -1L
    } else {
      x <- el[el[, 1] > n & el[, 2] > n, , drop = FALSE] - n
      if (length(x)) m[x] <- m[x[, 2:1, drop = FALSE]] <- -1L
      el <- el[el[, 1] <= n & el[, 2] <= n, , drop = FALSE]
    }
    if (length(el)) m[el] <- m[el[, 2:1, drop = FALSE]] <- 1L
  } else {
    m <- matrix(0L, n, n)
    if (length(el)) m[el] <- m[el[, 2:1, drop = FALSE]] <- 1L
  }
  m
}

#' Simulate from a fitted model
#'
#' Posterior predictive simulation of the formation and persistence networks
#' of every observed transition, given the observed interactions (sign
#' process) or the observed previous networks (interaction process), together
#' with the implied networks at each time point,
#' \eqn{y_t = y^P_t \cup (y^F_t \setminus y_{t-1})}.
#'
#' For the sign process, the simulated positive ties at time \eqn{t} are
#' \eqn{z_t = \{z^P_t \cup (z^F_t \setminus z_{t-1})\} \cap x_t}, and all
#' other interactions of \eqn{x_t} are negative.
#'
#' @param object a [bstergm()] or [bstergm_sign()] fit.
#' @param nsim number of simulations.
#' @param seed random seed.
#' @param coef optional parameter values (free parameters, as in
#'   `coef(object)`): a list with elements `F` and/or `P`. By default, `nsim`
#'   draws spread over the posterior sample are used.
#' @param output `"networks"` (signed or binary adjacency matrices) or
#'   `"stats"` (model statistics).
#' @param aux.iters number of Metropolis-Hastings steps of each simulation;
#'   by default five times the number used for fitting.
#' @param cores number of processes. Each simulation has its own random number
#'   stream, so results do not depend on the number of cores.
#' @param ... unused.
#'
#' @return For `output = "stats"`, a list with, for each component, an array
#'   of statistics (simulations x statistics x transitions). For `output =
#'   "networks"`, a list of `nsim` simulations, each a list with the simulated
#'   formation (`F`) and persistence (`P`) networks of every transition and,
#'   when both components were fitted, the implied networks `y` at times
#'   \eqn{1, \dots, T}.
#' @seealso [bgof()], [simulate_stats()]
#' @examples
#' data(senators70_am)
#' senate <- signed_series(senators70_am[1:3],
#'                         vertex.attr = attr(senators70_am$S99, "senators"))
#' fit <- bstergm_sign(~ edges + nodematch("Party") + gwesf_pos(0.6) + gwese_pos(0.6),
#'                data = senate, iters = 100, burn.in = 50, aux.iters = 1000, seed = 1,
#'                verbose = FALSE)
#' sims <- simulate(fit, nsim = 2, seed = 3)
#' sapply(sims[[1]]$y, signed_triad_census)
#' @export
simulate.bstergm <- function(object, nsim = 1, seed = NULL, coef = NULL,
                             output = c("networks", "stats"), aux.iters = NULL,
                             cores = getOption("Bstergm.cores", 1L), ...) {
  output <- match.arg(output)
  components <- object$components
  if (!is.null(coef)) {
    if (!is.list(coef)) {
      if (length(components) > 1) stop("'coef' must be a list with elements F and P.", call. = FALSE)
      coef <- stats::setNames(list(coef), components)
    }
    components <- intersect(components, names(coef))
  }
  Tn <- object$n_transitions
  n <- object$data$n
  formulas <- .formula_list(object$formula, Tn)
  ctrl <- object$control
  setup <- list()
  for (U in components) {
    ctrl$aux.iters <- aux.iters
    targets <- .build_targets(formulas, object$data, object$process, U, object$engine, ctrl)
    if (is.null(aux.iters)) for (t in seq_along(targets))
      targets[[t]]$aux.iters <- 5L * object$aux.iters[[U]][t]
    m <- as.ergm_model(targets[[1]]$state)
    draws <- if (is.null(coef)) .pick_draws(object, U, nsim)
             else matrix(rep(as.numeric(coef[[U]]), each = nsim), nsim)
    if (ncol(draws) != sum(object$free[[U]]))
      stop("'coef$", U, "' must have ", sum(object$free[[U]]), " values.", call. = FALSE)
    setup[[U]] <- list(targets = targets, etamap = m$etamap, draws = draws,
                       names = param_names(m, canonical = TRUE))
  }

  cores <- max(1L, as.integer(cores))
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  jobs <- .sim_jobs(components, nsim, cores, seed)

  run_job <- function(job) {
    st <- setup[[job$component]]
    lapply(seq_along(job$rows), function(i) {
      r <- job$rows[i]
      .use_stream(job$rngs[[i]])
      eta <- ergm.eta(.full_theta(object, job$component, st$draws[r, ]), st$etamap)
      sims <- .aux_stats(st$targets, eta, networks = output == "networks")
      if (output == "stats") return(vapply(sims, `[[`, numeric(length(st$names)), "stats"))
      lapply(seq_along(sims), function(t)
        .sim_matrix(sims[[t]]$el, st$targets[[t]], object$process, object$engine, n))
    })
  }
  res <- .run_jobs(jobs, run_job, cores)
  per <- lapply(stats::setNames(components, components), function(U)
    unlist(res[vapply(jobs, `[[`, "", "component") == U], recursive = FALSE))

  tl <- paste(head(object$data$times, -1), tail(object$data$times, -1), sep = "->")
  if (output == "stats") {
    return(lapply(stats::setNames(components, components), function(U) {
      arr <- array(NA_real_, c(nsim, length(setup[[U]]$names), Tn),
                   dimnames = list(NULL, setup[[U]]$names, tl))
      for (k in seq_len(nsim)) arr[k, , ] <- per[[U]][[k]]
      arr
    }))
  }

  dn <- list(object$data$node.names, object$data$node.names)
  prev <- object$data$y[seq_len(Tn)]
  cur <- object$data$y[seq_len(Tn) + 1]
  out <- lapply(seq_len(nsim), function(k) {
    sim <- lapply(stats::setNames(components, components), function(U) {
      l <- lapply(per[[U]][[k]], function(m) { dimnames(m) <- dn; m })
      stats::setNames(l, tl)
    })
    if (all(c("F", "P") %in% components)) {
      sim$y <- stats::setNames(lapply(seq_len(Tn), function(t) {
        if (object$process == "sign") {
          zF <- sim$F[[t]] > 0; zP <- sim$P[[t]] > 0
          x <- cur[[t]] != 0
          z <- (zP | (zF & !(prev[[t]] > 0))) & x
          .signed_from_layers(z, x, dn)
        } else {
          xF <- sim$F[[t]] > 0; xP <- sim$P[[t]] > 0
          m <- (xP | (xF & !(prev[[t]] != 0))) * 1L
          dimnames(m) <- dn
          m
        }
      }), object$data$times[-1])
    }
    sim
  })
  attr(out, "process") <- object$process
  out
}

#' Simulate the statistics of each transition
#'
#' Simulates the formation and persistence networks of every transition of a
#' sequence of signed networks for given parameter values and returns their
#' statistics. Each simulation starts at the observed formation or
#' persistence network.
#'
#' @param formula a one-sided ergm formula (see [bstergm()]).
#' @param data a [signed_series()] object or a list of signed networks.
#' @param coef_F,coef_P parameters of the formation and persistence models
#'   (excluding statistics that are constant on the sample space, such as the
#'   `L(~edges, ~active)` term of the multilayer formulation). Either may be
#'   `NULL` to skip that component.
#' @param process `"sign"` or `"interaction"`.
#' @param nsim number of simulations.
#' @param aux.iters number of Metropolis-Hastings steps of each simulation (see
#'   [control.bstergm()] for the default).
#' @param control a [control.bstergm()] list.
#' @param seed random seed.
#' @param ... arguments of version 0.1 of this function: `y_t` (used as
#'   `data`) and `ergm_sim_iters` (used as `aux.iters`) are translated;
#'   `L.x_F`, `L.x_P`, `L.z_F`, `L.z_P`, `times` and `model.dim` are no longer
#'   needed and are ignored, as are `network` and `network.stats` (use
#'   [simulate.bstergm()] for simulated networks).
#'
#' @return A list with elements `sz._F` and `sz._P`: matrices (statistics x
#'   transitions) for `nsim = 1`, or arrays (simulations x statistics x
#'   transitions) otherwise. Constant statistics are included.
#' @seealso [simulate.bstergm()]
#' @examples
#' data(senators70_am)
#' senate <- signed_series(senators70_am[1:3])
#' simulate_stats(~ edges + gwesf_pos(0.6) + gwese_pos(0.6), senate,
#'                coef_F = c(-2, 0.5, -0.2), coef_P = c(-2, 0.5, -0.2), seed = 1)
#' @export
simulate_stats <- function(formula, data, coef_F = NULL, coef_P = NULL, process = c("sign", "interaction"),
                           nsim = 1, aux.iters = NULL, control = control.bstergm(), seed = NULL, ...) {
  dots <- list(...)
  if (missing(data) && !is.null(dots$y_t)) data <- dots$y_t
  if (is.null(aux.iters) && !is.null(dots$ergm_sim_iters)) aux.iters <- dots$ergm_sim_iters
  if (isTRUE(dots$network) || isTRUE(dots$network.stats))
    warning("Arguments 'network' and 'network.stats' of version 0.1 are no longer used; ",
            "use simulate() on a bstergm fit to obtain simulated networks.", call. = FALSE)
  process <- match.arg(process)
  data <- signed_series(data)
  Tn <- length(data$y) - 1L
  formulas <- .formula_list(formula, Tn)
  engine <- if (process == "sign") .sign_engine(formulas) else .interaction_engine(formulas)
  if (!is.null(aux.iters)) control$aux.iters <- aux.iters
  coefs <- list(F = coef_F, P = coef_P)
  out <- list(sz._F = NULL, sz._P = NULL)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  streams <- .rng_streams(seed, 2)
  for (k in 1:2) {
    U <- c("F", "P")[k]
    if (is.null(coefs[[U]])) next
    targets <- .build_targets(formulas, data, process, U, engine, control)
    info <- .model_info(targets)
    free <- !info$offset
    if (!info$curved) free <- free & !.constant_stats(.mple_design(targets, info))
    if (length(coefs[[U]]) != sum(free))
      stop("'coef_", U, "' must have ", sum(free), " values (", paste(info$names[free], collapse = ", "), ").",
           call. = FALSE)
    th <- numeric(length(free))
    th[free] <- coefs[[U]]
    eta <- ergm.eta(th, info$etamap)
    sims <- .run_jobs(list(list(rng = streams[[k]])), function(job)
      lapply(seq_len(nsim), function(i) vapply(.aux_stats(targets, eta), `[[`,
                                                numeric(length(info$canonical)), "stats")))[[1]]
    arr <- array(unlist(sims), c(length(info$canonical), Tn, nsim),
                 dimnames = list(info$canonical, paste0("t", seq_len(Tn)), NULL))
    out[[paste0("sz._", U)]] <- if (nsim == 1) {
      matrix(arr[, , 1], nrow = dim(arr)[1], dimnames = dimnames(arr)[1:2])
    } else aperm(arr, c(3, 1, 2))
  }
  out
}
