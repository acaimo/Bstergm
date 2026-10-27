#' Tuning parameters for bstergm
#'
#' @param aux.iters number of Metropolis-Hastings steps of each auxiliary
#'   network simulation. If `NULL` (default) it is set separately for each
#'   transition to `aux.visits` times its number of free dyads (but at least
#'   `aux.min`), so that the quality of the auxiliary draws does not depend on
#'   network size. See [aux_diagnostics()] to check that it is large enough.
#' @param aux.visits,aux.min see `aux.iters`.
#' @param MCMC.prop,MCMC.prop.weights,MCMC.prop.args proposal of the auxiliary
#'   simulations, passed to [ergm::control.simulate.formula()]. The default
#'   `~.` toggles a uniformly chosen free dyad; `~sparse + .triadic` gives
#'   ergm's default triad-focused proposal.
#' @param translate if `TRUE` (default), sign-process formulas written with
#'   the multilayer terms of \pkg{ergm.multi} (`L()`, `gwespL()`), as in
#'   version 0.1, are translated to the equivalent native signed terms, which
#'   are much faster. The translation is used only if it reproduces the
#'   observed statistics exactly.
#' @param adapt.start,adapt.interval,adapt.stop adaptive Metropolis-Hastings:
#'   the proposal covariance is updated every `adapt.interval` iterations
#'   between iterations `adapt.start` and `adapt.stop`.
#' @param adapt.scale scaling of the empirical covariance; `NULL` uses
#'   \eqn{2.38^2/d}, where \eqn{d} is the number of parameters.
#' @param adapt.eps small constant added to the diagonal of the empirical
#'   covariance.
#' @param proposal.sd initial proposal standard deviation(s) of adaptive
#'   Metropolis-Hastings. If `NULL` (default) the initial proposal covariance
#'   is the scaled pseudo-likelihood covariance, or `0.025^2` times the
#'   identity when it is not available.
#' @param ads.gamma move factor of adaptive direction sampling.
#' @param ads.var variance of the Gaussian noise added to adaptive direction
#'   sampling proposals.
#' @param init.jitter scale of the random perturbation of the starting values:
#'   chains start at the pseudo-likelihood estimate plus `init.jitter` times a
#'   draw from its approximate covariance.
#' @param ads.init.iters,ads.init.chains length and population size (default
#'   `max(3, 2d)`) of the initial adaptive direction sampling run used with
#'   `init = "ads"` in [bstergm()] and [bstergm_sign()].
#' @param print.every report progress every so many iterations when `verbose`.
#'
#' @return A list of class `control.bstergm`.
#' @seealso [bstergm()], [bstergm_sign()]
#' @export
control.bstergm <- function(aux.iters = NULL, aux.visits = 3, aux.min = 500,
                            MCMC.prop = ~., MCMC.prop.weights = "default", MCMC.prop.args = list(),
                            translate = TRUE,
                            adapt.start = 100, adapt.interval = 50, adapt.stop = Inf,
                            adapt.scale = NULL, adapt.eps = 1e-6, proposal.sd = NULL,
                            ads.gamma = 0.5, ads.var = 0.025, init.jitter = 0.5,
                            ads.init.iters = 1000, ads.init.chains = NULL,
                            print.every = 1000) {
  if (inherits(MCMC.prop, "formula") && identical(environment(MCMC.prop), environment()))
    environment(MCMC.prop) <- globalenv()
  structure(as.list(environment()), class = "control.bstergm")
}

#' Bayesian separable temporal ERGMs
#'
#' `bstergm()` samples the posterior distribution of the formation and
#' persistence parameters of a separable temporal ERGM (STERGM) for a sequence
#' of binary networks. `bstergm_sign()` does the same for the sign process of
#' the separable two-layer STERGM of a sequence of signed networks. Both use
#' the approximate exchange algorithm with adaptive Metropolis-Hastings
#' (`method = "amh"`) or adaptive direction sampling (`method = "ads"`).
#'
#' @param formula a one-sided ergm formula, or a list of such formulas (one
#'   per transition, for time-varying covariates) with identical terms. For
#'   `bstergm()` it is an ordinary ergm formula evaluated on the formation and
#'   persistence networks. For `bstergm_sign()` it is evaluated on the
#'   positive ties of each formation or persistence network: ordinary ergm
#'   terms (e.g. `edges`, `nodematch("Party")`, `gwdegree(0.2, fixed = TRUE)`,
#'   `edgecov(x)`) refer to positive ties, and the signed-network terms of
#'   this package ([gwesf_pos][gwesf_pos-ergmTerm],
#'   [gwese_pos][gwese_pos-ergmTerm], [gwese_neg][gwese_neg-ergmTerm],
#'   [gwdegree_neg][gwdegree_neg-ergmTerm], ...) refer to positive and
#'   negative ties. The multilayer terms of \pkg{ergm.multi} on the layers
#'   `positive` and `active`, as in earlier versions of the package, are also
#'   accepted (slower).
#' @param data for `bstergm()`, a list of at least two undirected binary
#'   networks (adjacency matrices or `network` objects), a [network_series()]
#'   object, or a [signed_series()] object, whose interactions (ties of either
#'   sign) are then modelled. For `bstergm_sign()`, a [signed_series()]
#'   object, or a list of signed networks accepted by it.
#' @param method `"amh"` (adaptive Metropolis-Hastings, default) or `"ads"`
#'   (adaptive direction sampling).
#' @param components which components to fit: `"F"` (formation) and/or `"P"`
#'   (persistence).
#' @param iters total number of iterations of each chain, including burn-in.
#' @param burn.in number of initial iterations discarded.
#' @param thin thinning interval of the stored draws.
#' @param n_chains number of chains: independent chains for `"amh"` (default
#'   1), the population size for `"ads"` (default twice the number of
#'   parameters, at least 3).
#' @param aux.iters number of steps of each auxiliary network simulation; see
#'   [control.bstergm()] for the default.
#' @param prior.mean,prior.sigma mean and covariance of the multivariate
#'   normal prior of each component: scalars, vectors (diagonal covariance),
#'   matrices, or lists with elements `F` and `P` for component-specific
#'   priors. The defaults are \eqn{N(0, 10 I)}.
#' @param init starting values: `"mple"` (default) starts from the
#'   penalised maximum pseudo-likelihood estimate pooled over transitions;
#'   `"ads"` (with `method = "amh"`) first runs adaptive direction sampling
#'   for `control.bstergm(ads.init.iters =)` iterations, as in Caimo and
#'   Gollini, and starts the adaptive Metropolis-Hastings chains from its
#'   final population, with an initial proposal covariance estimated from it;
#'   `"prior"` starts from the prior mean plus a small uniform perturbation
#'   (the behaviour of version 0.1); a numeric vector, matrix (one column per
#'   chain), or list with elements `F` and `P`; or a previous fit of the same
#'   model, to continue its chains from their last stored states with the
#'   adapted proposal covariance.
#' @param offset.coef values of the parameters of `offset()` terms.
#' @param cores number of processes (default `getOption("Bstergm.cores",
#'   1)`). The formation and persistence components and, for `"amh"`, the
#'   chains are independent and are run in parallel (forking on Unix-alikes, a
#'   socket cluster on Windows). Results do not depend on the number of cores.
#' @param seed random seed, for reproducible results.
#' @param control a list of tuning parameters, see [control.bstergm()].
#' @param verbose print progress messages.
#'
#' @details
#' \strong{Binary networks.} For each transition \eqn{y_{t-1} \to y_t} the
#' formation network \eqn{y^F_t = y_{t-1} \cup y_t} and the persistence
#' network \eqn{y^P_t = y_{t-1} \cap y_t} are modelled by ERGMs with
#' parameters \eqn{\theta^F} on \eqn{\{y^F \supseteq y_{t-1}\}} and
#' \eqn{\theta^P} on \eqn{\{y^P \subseteq y_{t-1}\}}, as in the separable
#' temporal ERGM of Krivitsky and Handcock (2014).
#'
#' \strong{Signed networks.} A signed network is decomposed into an
#' interaction layer \eqn{x} and a sign layer \eqn{z} (see
#' [signed_series()]), and the formation and persistence networks are formed
#' for both layers. `bstergm_sign()` models \eqn{z^U_t \mid x^U_t}
#' (\eqn{U = F, P}) by an ERGM with parameter \eqn{\zeta^U} on the sample space
#' in which the sign of every interacting dyad of \eqn{x^U_t} may vary and the
#' interactions are fixed. The interaction layer follows a binary STERGM and
#' is fitted by `bstergm()` on the same [signed_series()]. Under independent
#' priors the posterior factorises over the components, which are sampled
#' independently (and in parallel).
#'
#' \strong{Algorithm.} The likelihood is intractable, so the posterior is
#' sampled with the approximate exchange algorithm (Caimo and Friel, 2011):
#' at each iteration a parameter \eqn{\theta'} is proposed, an auxiliary
#' network is simulated from the model at \eqn{\theta'} for every transition,
#' starting at the observed network, and \eqn{\theta'} is accepted with
#' probability
#' \deqn{\min\left\{1, \exp\left[(\eta(\theta) - \eta(\theta'))^\top \sum_t \{s(y'_t) - s(y_t)\}\right] \frac{\pi(\theta')}{\pi(\theta)}\right\}}{min{1, exp[(eta(theta) - eta(theta'))' sum_t {s(y'_t) - s(y_t)}] pi(theta') / pi(theta)}}
#' With `method = "amh"` proposals are Gaussian random walks whose covariance
#' is adapted to the chain history (Haario et al., 2001); with `method =
#' "ads"` each chain of a population moves along the difference of two other
#' chains (Caimo and Friel, 2011). Caimo and Gollini use adaptive direction
#' sampling only to find starting values during burn-in and then adaptive
#' Metropolis-Hastings; this scheme is obtained with `method = "amh"` and
#' `init = "ads"`.
#'
#' \strong{Computation.} The ergm model, proposal and constraints of every
#' transition are initialised once, and only the C sampler of \pkg{ergm} is
#' called during the iterations; with the signed-network terms of this
#' package the model statistics are also computed in C. Statistics that are
#' constant on the sample space (such as `L(~edges, ~active)` of the
#' multilayer formulation) are detected and dropped.
#'
#' \strong{Length of the auxiliary simulations.} The approximate exchange
#' algorithm treats the end of each auxiliary simulation as an exact draw from
#' the model. Simulations that are too short return networks close to the
#' observed ones, and the sampled distribution then resembles a
#' contrastive-divergence approximation more than the posterior. This matters
#' most for near-degenerate models. Check the auxiliary simulations with
#' [aux_diagnostics()] and the fit with [bgof()]. If needed, increase
#' `aux.iters` or `control.bstergm(aux.visits =)`; the cost per iteration
#' grows linearly.
#'
#' \strong{Lagged covariates.} The networks of each transition carry the
#' previous time point as network attributes, so that lagged dyadic
#' covariates can be used as `edgecov("lag.active")` (any tie at time
#' \eqn{t-1}) and, for signed data, `edgecov("lag.pos")` and
#' `edgecov("lag.neg")`.
#'
#' @return An object of class `bstergm`, a list including
#' \item{theta}{a list with, for each fitted component, an array of posterior
#'   draws (iterations x parameters x chains);}
#' \item{theta_F, theta_P}{the same arrays (for compatibility with version 0.1);}
#' \item{acceptance}{acceptance rates of each chain;}
#' \item{mple}{the pseudo-likelihood estimates used for initialisation;}
#' \item{init.ads}{with `init = "ads"`, a summary of the initial adaptive
#'   direction sampling run of each component;}
#' \item{obs.stats}{the observed statistics of each transition;}
#' \item{process}{`"sign"` for `bstergm_sign()`, `"interaction"` for
#'   `bstergm()`;}
#' \item{prior, control, formula, data, time}{inputs and timings.}
#' Use [summary.bstergm()], [plot.bstergm()], [bgof()] and
#' [as.mcmc.list.bstergm()] to analyse it.
#'
#' @references
#' Caimo, A. and Gollini, I. (2025). Separable models for dynamic signed
#' networks. arXiv:2505.07669.
#'
#' Caimo, A. and Friel, N. (2011). Bayesian inference for exponential random
#' graph models. \emph{Social Networks}, 33(1), 41-55.
#'
#' Haario, H., Saksman, E. and Tamminen, J. (2001). An adaptive Metropolis
#' algorithm. \emph{Bernoulli}, 7(2), 223-242.
#'
#' Krivitsky, P. N. and Handcock, M. S. (2014). A separable model for dynamic
#' networks. \emph{Journal of the Royal Statistical Society B}, 76(1), 29-46.
#'
#' @seealso [network_series()], [signed_series()], [summary.bstergm()],
#'   [bgof()], [aux_diagnostics()], [control.bstergm()], [predict.bstergm()]
#' @examples
#' data(senators70_am)
#' senate <- signed_series(senators70_am[1:3],
#'                         vertex.attr = attr(senators70_am$S99, "senators"))
#'
#' # Very short runs for illustration; in practice use many more iterations
#' # (e.g. iters = 30000, burn.in = 10000) and the default auxiliary length.
#'
#' # Binary STERGM of the interactions (ties of either sign)
#' fit_x <- bstergm(~ edges + nodematch("Party"), data = senate,
#'                  iters = 100, burn.in = 50, aux.iters = 1000, seed = 1, verbose = FALSE)
#' summary(fit_x)
#'
#' # Sign process given the interactions
#' fit_z <- bstergm_sign(~ edges + nodematch("Party") + gwesf_pos(0.6) + gwese_pos(0.6),
#'                       data = senate, iters = 100, burn.in = 50, aux.iters = 1000,
#'                       seed = 1, verbose = FALSE)
#' summary(fit_z)
#' @export
bstergm <- function(formula, data, method = c("amh", "ads"), components = c("F", "P"),
                    iters = 30000, burn.in = 10000, thin = 1, n_chains = NULL,
                    aux.iters = NULL, prior.mean = 0, prior.sigma = 10,
                    init = "mple", offset.coef = NULL,
                    cores = getOption("Bstergm.cores", 1L), seed = NULL,
                    control = control.bstergm(), verbose = TRUE) {
  args <- as.list(environment())
  args$data <- network_series(data)
  out <- do.call(.bstergm_fit, c(args, list(process = "interaction")))
  out$call <- match.call()
  out
}

#' @rdname bstergm
#' @export
bstergm_sign <- function(formula, data, method = c("amh", "ads"), components = c("F", "P"),
                         iters = 30000, burn.in = 10000, thin = 1, n_chains = NULL,
                         aux.iters = NULL, prior.mean = 0, prior.sigma = 10,
                         init = "mple", offset.coef = NULL,
                         cores = getOption("Bstergm.cores", 1L), seed = NULL,
                         control = control.bstergm(), verbose = TRUE) {
  args <- as.list(environment())
  out <- do.call(.bstergm_fit, c(args, list(process = "sign")))
  out$call <- match.call()
  out
}

.bstergm_fit <- function(formula, data, process, method = c("amh", "ads"), components = c("F", "P"),
                         iters = 30000, burn.in = 10000, thin = 1, n_chains = NULL,
                         aux.iters = NULL, prior.mean = 0, prior.sigma = 10,
                         init = "mple", offset.coef = NULL,
                         cores = getOption("Bstergm.cores", 1L), seed = NULL,
                         control = control.bstergm(), verbose = TRUE) {
  cl <- NULL
  t.start <- proc.time()[[3]]
  method <- match.arg(method)
  components <- unique(match.arg(components, c("F", "P"), several.ok = TRUE))
  if (!inherits(control, "control.bstergm")) control <- do.call(control.bstergm, as.list(control))
  if (!is.null(aux.iters)) control$aux.iters <- aux.iters
  iters <- as.integer(iters); burn.in <- as.integer(burn.in); thin <- max(1L, as.integer(thin))
  if (iters < 1 || burn.in < 0) stop("'iters' must be positive and 'burn.in' non-negative.", call. = FALSE)
  if (burn.in >= iters) stop("'burn.in' must be smaller than 'iters'.", call. = FALSE)

  data <- signed_series(data)
  Tn <- length(data$y) - 1L
  formulas <- .formula_list(formula, Tn)
  engine <- if (process == "sign") .sign_engine(formulas) else .interaction_engine(formulas)
  formula.original <- formula
  if (engine == "layer" && isTRUE(control$translate)) {
    tr <- .translate_formulas(formulas, data)
    if (!is.null(tr)) {
      formulas <- tr
      formula <- if (length(unique(lapply(tr, deparse))) == 1 && inherits(formula.original, "formula")) tr[[1]] else tr
      engine <- "native"
      if (verbose) message("Using the equivalent native signed terms: ",
                           paste(deparse(tr[[1]], width.cutoff = 500L), collapse = " "))
    }
  }

  if (verbose) message("Initialising the ", process, " process models (", engine, " engine, ",
                       Tn, " transition", if (Tn > 1) "s", ") ...")
  comp <- list()
  for (U in components) {
    targets <- .build_targets(formulas, data, process, U, engine, control)
    comp[[U]] <- .setup_component(targets, U, method, prior.mean, prior.sigma, init,
                                  offset.coef, n_chains, control, verbose)
  }

  jobs <- list()
  for (U in components) {
    H <- comp[[U]]$n_chains
    if (method == "amh") {
      for (h in seq_len(H)) jobs[[length(jobs) + 1]] <- list(component = U, chains = h)
    } else jobs[[length(jobs) + 1]] <- list(component = U, chains = seq_len(H))
  }
  pre <- Filter(function(U) isTRUE(comp[[U]]$ads.init), components)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  streams <- .rng_streams(seed, length(jobs) + length(pre))
  for (j in seq_along(jobs)) jobs[[j]]$rng <- streams[[j]]

  ctl <- c(control, list(iters = iters, burn.in = burn.in, thin = thin))
  nk <- length(.keep_iterations(iters, burn.in, thin))

  if (length(pre)) {
    pjobs <- lapply(seq_along(pre), function(k) list(component = pre[k], rng = streams[[length(jobs) + k]]))
    if (verbose) message("Initial adaptive direction sampling (", control$ads.init.iters,
                         " iterations) for component", if (length(pre) > 1) "s", " ",
                         paste(pre, collapse = " and "), " ...")
    pres <- .run_jobs(pjobs, function(job) {
      U <- job$component
      cp <- comp[[U]]
      H0 <- NVL(control$ads.init.chains, max(3L, 2L * cp$spec$p))
      ctl0 <- ctl
      ctl0$iters <- max(2L, as.integer(control$ads.init.iters))
      ctl0$burn.in <- ctl0$iters %/% 2L
      ctl0$thin <- 1L
      ctl0$ads.gamma <- .pick(ctl$ads.gamma, U)
      ctl0$ads.var <- .pick(ctl$ads.var, U)
      .run_ads(cp$targets, cp$spec, ctl0, .init_values(cp, seq_len(H0), control),
               paste(U, "initial ADS"), verbose)
    }, cores)
    for (k in seq_along(pre)) {
      comp[[pre[k]]] <- .ads_start(comp[[pre[k]]], pres[[k]], control)
    }
  }
  if (verbose) {
    message("Sampling: ", paste(components, collapse = " and "), " component",
            if (length(components) > 1) "s", ", ", method, ", ",
            iters, " iterations (", burn.in, " burn-in), ",
            length(jobs), " job", if (length(jobs) > 1) "s", " on ",
            max(1L, min(cores, length(jobs))), " core", if (min(cores, length(jobs)) > 1) "s", ".")
    for (U in components)
      message("  ", U, ": auxiliary simulations of ",
              paste(vapply(comp[[U]]$targets, `[[`, 1L, "aux.iters"), collapse = ", "), " steps")
  }

  run_job <- function(job) {
    cp <- comp[[job$component]]
    lab <- if (method == "amh") sprintf("%s chain %d", job$component, job$chains)
           else sprintf("%s population", job$component)
    init.mat <- .init_values(cp, job$chains, control)
    ctlU <- ctl
    ctlU$ads.gamma <- .pick(ctl$ads.gamma, job$component)
    ctlU$ads.var <- .pick(ctl$ads.var, job$component)
    if (method == "amh")
      .run_amh(cp$targets, cp$spec, ctlU, init.mat[, 1], cp$prop.cov, lab, verbose)
    else
      .run_ads(cp$targets, cp$spec, ctlU, init.mat, lab, verbose)
  }
  res <- .run_jobs(jobs, run_job, cores)

  theta <- acceptance <- proposal <- jtime <- list()
  for (U in components) {
    sp <- comp[[U]]$spec
    idx <- which(vapply(jobs, `[[`, "", "component") == U)
    if (method == "amh") {
      H <- length(idx)
      arr <- array(NA_real_, c(nk, sp$p, H), dimnames = list(NULL, sp$names, paste0("chain", seq_len(H))))
      for (h in seq_len(H)) arr[, , h] <- res[[idx[h]]]$theta
      theta[[U]] <- arr
      acceptance[[U]] <- vapply(res[idx], `[[`, 0, "accept")
      proposal[[U]] <- lapply(res[idx], `[[`, "proposal")
    } else {
      theta[[U]] <- res[[idx]]$theta
      acceptance[[U]] <- res[[idx]]$accept
    }
    names(acceptance[[U]]) <- paste0("chain", seq_along(acceptance[[U]]))
    jtime[[U]] <- vapply(res[idx], `[[`, 0, "time")
  }

  out <- list(
    call = cl, formula = formula, formula.original = formula.original,
    process = process, method = method, engine = engine,
    components = components, theta = theta,
    theta_F = theta$F, theta_P = theta$P,
    acceptance = acceptance, proposal = proposal,
    coef.names = lapply(comp, function(cp) cp$spec$names),
    fixed = lapply(comp, function(cp) cp$fixed),
    free = lapply(comp, function(cp) cp$spec$free),
    theta.fixed = lapply(comp, function(cp) stats::setNames(cp$spec$fixed, cp$info$names)),
    mple = lapply(comp, `[[`, "mple"),
    init.ads = Filter(Negate(is.null), lapply(comp, `[[`, "ads.run")),
    prior = lapply(comp, function(cp) list(mean = cp$spec$prior$mean, sigma = cp$prior.sigma)),
    obs.stats = lapply(comp, function(cp) .obs_matrix(cp$targets, cp$info$canonical)),
    aux.iters = lapply(comp, function(cp) vapply(cp$targets, `[[`, 1L, "aux.iters")),
    n_chains = vapply(comp, `[[`, 1L, "n_chains"),
    iters = iters, burn.in = burn.in, thin = thin, seed = seed,
    control = control, data = data, n_transitions = Tn,
    time = list(total = proc.time()[[3]] - t.start, sampling = jtime))
  class(out) <- "bstergm"
  if (verbose) message("Done in ", .format_time(out$time$total), ".")
  out
}

.interaction_engine <- function(formulas) {
  nms <- unique(unlist(lapply(formulas, all.names)))
  bad <- intersect(nms, c(.signed_term_names, .multilayer_term_names()))
  if (length(bad))
    stop("bstergm() uses ordinary ergm terms; ", paste(sQuote(bad), collapse = ", "),
         " need bstergm_sign().", call. = FALSE)
  "standard"
}

.obs_matrix <- function(targets, names) {
  m <- do.call(rbind, lapply(targets, `[[`, "obs"))
  dimnames(m) <- list(paste0("t", seq_len(nrow(m))), names)
  m
}

.pick <- function(x, U) if (is.list(x) && !is.null(x[[U]])) x[[U]] else x

.prior_component <- function(mean, sigma, U, p, names, free) {
  m <- .pick(mean, U); S <- .pick(sigma, U)
  ptot <- length(free)
  if (is.null(m)) m <- 0
  if (is.null(S)) S <- 10
  m <- as.numeric(m)
  if (length(m) == ptot && ptot != p) m <- m[free]
  if (length(m) == 1) m <- rep(m, p)
  if (length(m) != p) stop("The prior mean of component ", U, " must have length ", p, ".", call. = FALSE)
  if (is.matrix(S)) {
    if (nrow(S) == ptot && ptot != p) S <- S[free, free, drop = FALSE]
  } else {
    S <- as.numeric(S)
    if (length(S) == ptot && ptot != p) S <- S[free]
    if (length(S) == 1) S <- rep(S, p)
    S <- diag(S, p)
  }
  if (!all(dim(S) == p)) stop("The prior covariance of component ", U, " must be ", p, " x ", p, ".", call. = FALSE)
  prec <- tryCatch(solve(S), error = function(e) stop("The prior covariance of component ", U,
                                                     " is not invertible.", call. = FALSE))
  names(m) <- names
  dimnames(S) <- list(names, names)
  list(mean = m, sigma = S, prec = prec)
}

.setup_component <- function(targets, U, method, prior.mean, prior.sigma, init, offset.coef,
                             n_chains, control, verbose) {
  info <- .model_info(targets)
  ptot <- length(info$names)
  free <- !info$offset
  fixed <- numeric(ptot)
  if (any(info$offset)) {
    if (is.null(offset.coef) || length(offset.coef) != sum(info$offset))
      stop("'offset.coef' must give a value for each of the ", sum(info$offset), " offset parameter(s).",
           call. = FALSE)
    if (any(!is.finite(offset.coef)))
      stop("Infinite offsets are not supported; use constraints instead.", call. = FALSE)
    fixed[info$offset] <- offset.coef
  }

  design <- NULL
  if (!info$curved) {
    design <- .mple_design(targets, info)
    const <- .constant_stats(design) & free
    if (any(const)) {
      if (verbose) message("  ", U, ": dropping statistic(s) constant on the sample space: ",
                           paste(info$names[const], collapse = ", "))
      free[const] <- FALSE
    }
  }
  p <- sum(free)
  if (p == 0) stop("Component ", U, " has no free parameters.", call. = FALSE)
  names.free <- info$names[free]
  prior <- .prior_component(prior.mean, prior.sigma, U, p, names.free, free)

  H <- if (!is.null(n_chains)) as.integer(n_chains) else if (method == "amh") 1L else max(3L, 2L * p)
  if (method == "ads" && H < 3) stop("Adaptive direction sampling needs 'n_chains' >= 3.", call. = FALSE)

  mple <- NULL
  prev.prop <- NULL
  if (inherits(init, "bstergm")) {
    if (is.null(init$theta[[U]]) || !identical(dimnames(init$theta[[U]])[[2]], names.free))
      stop("The previous fit given as 'init' does not match component ", U, ".", call. = FALSE)
    a <- init$theta[[U]]
    init.U <- matrix(a[dim(a)[1], , ], nrow = p)
    if (!is.null(init$proposal[[U]])) prev.prop <- init$proposal[[U]][[1]]
  } else init.U <- .pick(init, U)
  ads.init <- is.character(init.U) && identical(init.U, "ads")
  if (ads.init) {
    if (method != "amh") stop("'init = \"ads\"' is only available with method = \"amh\".", call. = FALSE)
    init.U <- "mple"
  }
  if (is.character(init.U) && identical(init.U, "mple")) {
    if (!is.null(design)) mple <- .mple_init(design, free, fixed, prior)
    if (is.null(mple) && verbose)
      message("  ", U, ": pseudo-likelihood initialisation not available; starting from the prior mean.")
  } else if (is.numeric(init.U)) {
    im <- as.matrix(init.U)
    if (nrow(im) == ptot && ptot != p) im <- im[free, , drop = FALSE]
    else if (ncol(im) == ptot && nrow(im) != p && ptot != p) im <- t(im[, free, drop = FALSE])
    else if (nrow(im) != p && ncol(im) == p) im <- t(im)
    if (nrow(im) != p || anyNA(im))
      stop("Numeric starting values of component ", U, " must have ", p, " values per chain (",
           paste(names.free, collapse = ", "), ").", call. = FALSE)
    init.U <- im
  } else if (!(is.character(init.U) && identical(init.U, "prior"))) {
    stop("'init' must be \"mple\", \"ads\", \"prior\", a numeric vector or matrix, a previous fit, ",
         "or a list of these.", call. = FALSE)
  }

  psd <- .pick(control$proposal.sd, U)
  prop.cov <- if (!is.null(psd)) diag(rep_len(psd, p)^2, p)
              else if (!is.null(prev.prop)) prev.prop
              else if (!is.null(mple)) (2.38^2 / p) * mple$cov
              else diag(0.025^2, p)

  spec <- list(p = p, names = names.free, free = free, fixed = fixed,
               etamap = info$etamap, prior = prior[c("mean", "prec")])
  list(targets = targets, info = info, spec = spec, prior.sigma = prior$sigma,
       fixed = if (any(!free)) setNames(fixed[!free], info$names[!free]) else NULL,
       mple = mple, init = init.U, n_chains = H, prop.cov = prop.cov, U = U,
       ads.init = ads.init, user.prop = !is.null(psd))
}

.ads_start <- function(cp, run, control) {
  p <- cp$spec$p
  last <- run$last
  H0 <- ncol(last)
  cp$init <- last[, rep_len(seq_len(H0), cp$n_chains), drop = FALSE]
  rownames(cp$init) <- cp$spec$names
  draws <- matrix(aperm(run$theta, c(1, 3, 2)), ncol = p)
  S <- if (nrow(draws) > p) stats::cov(draws) else NULL
  if (!cp$user.prop && !is.null(S)) {
    C <- NVL(control$adapt.scale, 2.38^2 / p) * (S + control$adapt.eps * diag(p))
    if (!is.null(.safe_chol(C))) cp$prop.cov <- C
  }
  if (!is.null(S)) dimnames(S) <- list(cp$spec$names, cp$spec$names)
  cp$ads.run <- list(iters = max(2L, as.integer(control$ads.init.iters)), chains = H0,
                     accept = mean(c(run$accept, run$accept.burnin)),
                     mean = colMeans(draws), cov = S, last = last, time = run$time)
  cp
}

.init_values <- function(cp, chains, control) {
  p <- cp$spec$p
  H <- length(chains)
  init <- cp$init
  if (is.numeric(init)) {
    return(init[, rep_len(seq_len(ncol(init)), max(chains))[chains], drop = FALSE])
  }
  if (!is.null(cp$mple)) {
    R <- chol(cp$mple$cov)
    z <- matrix(stats::rnorm(p * H), H, p) %*% R
    return(cp$mple$coef + control$init.jitter * t(z))
  }
  matrix(cp$spec$prior$mean + stats::runif(p * H, -0.1, 0.1), p, H)
}
