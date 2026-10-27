.eta <- function(spec, theta) {
  th <- spec$fixed
  th[spec$free] <- theta
  ergm.eta(th, spec$etamap)
}

.logprior <- function(spec, theta) {
  d <- theta - spec$prior$mean
  -0.5 * sum(d * drop(spec$prior$prec %*% d))
}

.safe_chol <- function(S) {
  tryCatch(chol(S), error = function(e) NULL)
}

.keep_iterations <- function(iters, burn.in, thin) {
  if (iters <= burn.in) return(integer(0))
  seq.int(burn.in + thin, iters, by = thin)
}

.progress <- function(label, i, iters, acc, nacc, t0) {
  el <- proc.time()[[3]] - t0
  message(sprintf("%s: iteration %d/%d, acceptance rate %.3f, %s elapsed, about %s left", label, i, iters,
                  acc / max(nacc, 1), .format_time(el), .format_time(el * (iters - i) / i)))
}

.format_time <- function(s) {
  if (s < 60) sprintf("%.1fs", s)
  else if (s < 3600) sprintf("%.1fmin", s / 60)
  else sprintf("%.2fh", s / 3600)
}

.run_amh <- function(targets, spec, ctl, init, prop.cov, label = "chain", verbose = FALSE) {
  p <- spec$p
  keep <- .keep_iterations(ctl$iters, ctl$burn.in, ctl$thin)
  out <- matrix(NA_real_, length(keep), p, dimnames = list(NULL, spec$names))
  theta <- init
  eta <- .eta(spec, theta)
  lp <- .logprior(spec, theta)
  C <- prop.cov
  R <- .safe_chol(C)
  if (is.null(R)) stop("The initial proposal covariance is not positive definite.", call. = FALSE)
  sd.scale <- NVL(ctl$adapt.scale, 2.38^2 / p)
  mu <- numeric(p); M2 <- matrix(0, p, p); nh <- 0
  acc <- 0L; acc.burn <- 0L; k <- 0L
  t0 <- proc.time()[[3]]
  for (i in seq_len(ctl$iters)) {
    nh <- nh + 1
    dl <- theta - mu
    mu <- mu + dl / nh
    M2 <- M2 + tcrossprod(dl, theta - mu)
    if (i > ctl$adapt.start && i <= ctl$adapt.stop && i %% ctl$adapt.interval == 0 && nh > 1) {
      Cn <- sd.scale * (M2 / (nh - 1) + ctl$adapt.eps * diag(p))
      Rn <- .safe_chol(Cn)
      if (!is.null(Rn)) { C <- Cn; R <- Rn }
    }
    prop <- theta + drop(stats::rnorm(p) %*% R)
    eta.p <- .eta(spec, prop)
    lp.p <- .logprior(spec, prop)
    delta <- .aux_delta(targets, eta.p)
    lr <- sum((eta - eta.p) * delta) + lp.p - lp
    if (is.finite(lr) && log(stats::runif(1)) < lr) {
      theta <- prop; eta <- eta.p; lp <- lp.p
      if (i > ctl$burn.in) acc <- acc + 1L else acc.burn <- acc.burn + 1L
    }
    if (k < length(keep) && i == keep[k + 1L]) {
      k <- k + 1L
      out[k, ] <- theta
    }
    if (verbose && i %% ctl$print.every == 0) {
      if (i > ctl$burn.in) .progress(label, i, ctl$iters, acc, i - ctl$burn.in, t0)
      else .progress(paste(label, "(burn-in)"), i, ctl$iters, acc.burn, i, t0)
    }
  }
  list(theta = out, accept = acc / max(ctl$iters - ctl$burn.in, 1),
       accept.burnin = acc.burn / max(min(ctl$burn.in, ctl$iters), 1),
       proposal = C, last = theta, time = proc.time()[[3]] - t0)
}

.run_ads <- function(targets, spec, ctl, init, label = "ADS", verbose = FALSE) {
  p <- spec$p
  H <- ncol(init)
  if (H < 3) stop("Adaptive direction sampling needs at least 3 chains.", call. = FALSE)
  keep <- .keep_iterations(ctl$iters, ctl$burn.in, ctl$thin)
  out <- array(NA_real_, c(length(keep), p, H),
               dimnames = list(NULL, spec$names, paste0("chain", seq_len(H))))
  theta <- init
  eta <- lapply(seq_len(H), function(h) .eta(spec, theta[, h]))
  lp <- vapply(seq_len(H), function(h) .logprior(spec, theta[, h]), 0)
  sdn <- sqrt(ctl$ads.var)
  acc <- integer(H); acc.burn <- integer(H); k <- 0L
  t0 <- proc.time()[[3]]
  for (i in seq_len(ctl$iters)) {
    for (h in seq_len(H)) {
      o <- seq_len(H)[-h][sample.int(H - 1L, 2L)]
      prop <- theta[, h] + ctl$ads.gamma * (theta[, o[2]] - theta[, o[1]]) + stats::rnorm(p, 0, sdn)
      eta.p <- .eta(spec, prop)
      lp.p <- .logprior(spec, prop)
      delta <- .aux_delta(targets, eta.p)
      lr <- sum((eta[[h]] - eta.p) * delta) + lp.p - lp[h]
      if (is.finite(lr) && log(stats::runif(1)) < lr) {
        theta[, h] <- prop; eta[[h]] <- eta.p; lp[h] <- lp.p
        if (i > ctl$burn.in) acc[h] <- acc[h] + 1L else acc.burn[h] <- acc.burn[h] + 1L
      }
    }
    if (k < length(keep) && i == keep[k + 1L]) {
      k <- k + 1L
      out[k, , ] <- theta
    }
    if (verbose && i %% ctl$print.every == 0) {
      if (i > ctl$burn.in) .progress(label, i, ctl$iters, sum(acc), H * (i - ctl$burn.in), t0)
      else .progress(paste(label, "(burn-in)"), i, ctl$iters, sum(acc.burn), H * i, t0)
    }
  }
  list(theta = out, accept = acc / max(ctl$iters - ctl$burn.in, 1),
       accept.burnin = acc.burn / max(min(ctl$burn.in, ctl$iters), 1),
       last = theta, time = proc.time()[[3]] - t0)
}
