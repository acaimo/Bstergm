.multilayer_term_names <- function() {
  ns <- asNamespace("ergm.multi")
  nm <- sub("^InitErgmTerm\\.", "", ls(ns, all.names = TRUE, pattern = "^InitErgmTerm\\."))
  unique(c("L", "CMBL", nm[grepl("L$", nm)]))
}

.signed_term_names <- c("gwesf_pos", "gwesf_neg", "gwese_pos", "gwese_neg", "gwesp_signed",
                        "esp_signed", "gwdegree_neg", "degree_neg", "triangle_signed")

.sign_engine <- function(formulas) {
  nms <- unique(unlist(lapply(formulas, all.names)))
  ml <- intersect(nms, .multilayer_term_names())
  if (length(ml) && length(intersect(nms, .signed_term_names)))
    stop("A formula cannot mix the multilayer terms of ergm.multi (", paste(ml, collapse = ", "),
         ") with the signed-network terms of Bstergm.", call. = FALSE)
  if (length(ml)) "layer" else "native"
}

.formula_list <- function(formula, Tn) {
  fl <- if (inherits(formula, "formula")) list(formula) else formula
  if (!is.list(fl) || !all(vapply(fl, inherits, TRUE, "formula")))
    stop("'formula' must be a formula or a list of formulas.", call. = FALSE)
  if (length(fl) == 1) fl <- rep(fl, Tn)
  if (length(fl) != Tn)
    stop("A list of formulas must have one formula per transition (", Tn, ").", call. = FALSE)
  lapply(fl, function(f) {
    if (length(f) == 3) f <- f[-2]
    f
  })
}

.binary_network <- function(m, vertex.attr = NULL) {
  nw <- .network_from_el(which(upper.tri(m) & m != 0, arr.ind = TRUE), nrow(m))
  if (!is.null(rownames(m))) network.vertex.names(nw) <- rownames(m)
  .set_vertex_attr(nw, vertex.attr)
}

.add_lags <- function(nw, prev) {
  nw %n% "lag.pos" <- .binary_network(prev > 0)
  nw %n% "lag.neg" <- .binary_network(prev < 0)
  nw %n% "lag.active" <- .binary_network(prev != 0)
  nw
}

.target_spec <- function(tn, t, process, component, engine, vertex.attr) {
  m <- tn[[component]][[t]]
  prev <- tn$prev[[t]]
  n <- nrow(m)
  if (process == "sign") {
    nfree <- sum(m != 0) / 2
    if (engine == "native") {
      basis <- .add_lags(.signed_network_from_matrix(m, vertex.attr), prev)
      constraints <- ~active_signs
    } else {
      z <- .binary_network(m > 0, vertex.attr)
      x <- .binary_network(m != 0, vertex.attr)
      basis <- Layer(positive = z, active = x)
      empty <- network.initialize(n, directed = FALSE)
      network.vertex.names(empty) <- network.vertex.names(x)
      constraints <- ~Dyads(vary = ~dyadcov(Layer(x, empty)))
    }
  } else {
    basis <- .add_lags(.binary_network(m != 0, vertex.attr), prev)
    prevx <- .binary_network(prev != 0)
    if (component == "F") {
      nfree <- n * (n - 1) / 2 - sum(prev != 0) / 2
      constraints <- if (network.edgecount(prevx)) ~fixedas(present = prevx) else ~.
    } else {
      nfree <- sum(prev != 0) / 2
      constraints <- ~fixallbut(free.dyads = prevx)
    }
  }
  list(basis = basis, constraints = constraints, nfree = nfree)
}

.aux_iters <- function(control, nfree) {
  if (!is.null(control$aux.iters)) return(as.integer(control$aux.iters))
  as.integer(max(control$aux.min, ceiling(control$aux.visits * nfree)))
}

.mcmc_control <- function(control, aux.iters) {
  control.simulate.formula(MCMC.burnin = aux.iters, MCMC.interval = 1L,
                           MCMC.prop = control$MCMC.prop,
                           MCMC.prop.weights = control$MCMC.prop.weights,
                           MCMC.prop.args = control$MCMC.prop.args)
}

.make_target <- function(formula, spec, control, monitor = NULL) {
  f <- if (is.null(monitor)) formula
       else statnet.common::append_rhs.formula(formula, list_rhs.formula(monitor), keep.onesided = TRUE)
  aux <- .aux_iters(control, spec$nfree)
  ctrl <- .mcmc_control(control, aux)
  a <- suppressMessages(
    simulate(f, coef = NULL, constraints = spec$constraints, basis = spec$basis,
             control = ctrl, output = "stats", return.args = "ergm_state"))
  state <- a$object
  list(state = ergm_state_send(state), control = a$control, obs = as.numeric(state$stats),
       aux.iters = aux, nfree = spec$nfree, formula = f, basis = spec$basis,
       constraints = spec$constraints)
}

.build_targets <- function(formulas, data, process, component, engine, control, monitor = NULL) {
  tn <- transition_networks(data)
  lapply(seq_along(formulas), function(t) {
    spec <- .target_spec(tn, t, process, component, engine, data$vertex.attr)
    .make_target(formulas[[t]], spec, control, monitor)
  })
}

.model_info <- function(targets) {
  m <- as.ergm_model(targets[[1]]$state)
  nms <- param_names(m, canonical = FALSE)
  cnms <- param_names(m, canonical = TRUE)
  for (tg in targets[-1]) {
    m2 <- as.ergm_model(tg$state)
    if (!identical(param_names(m2, canonical = FALSE), nms))
      stop("All transitions must have the same model terms.", call. = FALSE)
  }
  list(etamap = m$etamap, names = nms, canonical = cnms, curved = is.curved(m),
       offset = m$etamap$offsettheta)
}

.aux_delta <- function(targets, eta) {
  d <- 0
  for (tg in targets) {
    z <- ergm_MCMC_slave(tg$state, eta = eta, control = tg$control, verbose = FALSE,
                         burnin = tg$aux.iters, samplesize = 1L, interval = 1L)
    if (!is.null(z$status) && z$status != 0)
      stop("The auxiliary network simulation failed (ergm status ", z$status, ").", call. = FALSE)
    d <- d + (z$s[1, ] - tg$obs)
  }
  d
}

.aux_stats <- function(targets, eta, burnin = NULL, networks = FALSE) {
  lapply(targets, function(tg) {
    ctrl <- tg$control
    if (networks) ctrl$MCMC.save_networks <- TRUE
    z <- ergm_MCMC_slave(tg$state, eta = eta, control = ctrl, verbose = FALSE,
                         burnin = NVL(burnin, tg$aux.iters), samplesize = 1L, interval = 1L)
    if (!is.null(z$status) && z$status != 0)
      stop("The auxiliary network simulation failed (ergm status ", z$status, ").", call. = FALSE)
    out <- list(stats = z$s[1, ])
    if (networks) out$el <- z$state$el
    out
  })
}
