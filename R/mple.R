.mple_design <- function(targets, info) {
  parts <- lapply(targets, function(tg) {
    ergmMPLE(tg$formula, basis = tg$basis, constraints = tg$constraints,
             output = "matrix")
  })
  X <- do.call(rbind, lapply(parts, `[[`, "predictor"))
  y <- unlist(lapply(parts, `[[`, "response"))
  w <- unlist(lapply(parts, `[[`, "weights"))
  o <- do.call(order, c(unname(as.data.frame(X)), list(y, w)))
  X <- X[o, , drop = FALSE]
  colnames(X) <- info$canonical
  list(X = X, y = y[o], w = w[o])
}

.constant_stats <- function(design) {
  apply(design$X, 2, function(v) all(v == 0))
}

.penalised_logit <- function(X, y, w, off, m, Q, maxit = 100, tol = 1e-8) {
  b <- m
  for (it in seq_len(maxit)) {
    eta <- drop(X %*% b) + off
    mu <- stats::plogis(eta)
    g <- crossprod(X, w * (y - mu)) - Q %*% (b - m)
    H <- crossprod(X, X * (w * mu * (1 - mu))) + Q
    step <- tryCatch(solve(H, g), error = function(e) NULL)
    if (is.null(step)) break
    b <- b + drop(step)
    if (max(abs(step)) < tol) break
  }
  eta <- drop(X %*% b) + off
  mu <- stats::plogis(eta)
  H <- crossprod(X, X * (w * mu * (1 - mu))) + Q
  V <- tryCatch(solve(H), error = function(e) NULL)
  list(coef = b, cov = V, iter = it, converged = it < maxit)
}

.mple_init <- function(design, free, fixed.values, prior) {
  Xf <- design$X[, free, drop = FALSE]
  off <- if (any(!free)) drop(design$X[, !free, drop = FALSE] %*% fixed.values[!free]) else 0
  fit <- .penalised_logit(Xf, design$y, design$w, off, prior$mean, prior$prec)
  if (is.null(fit$cov) || any(!is.finite(fit$coef))) return(NULL)
  V <- (fit$cov + t(fit$cov)) / 2
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev <= 0)) return(NULL)
  names(fit$coef) <- colnames(Xf)
  dimnames(V) <- list(colnames(Xf), colnames(Xf))
  list(coef = fit$coef, cov = V, converged = fit$converged)
}
