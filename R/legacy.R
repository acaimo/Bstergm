.legacy_args <- function(args, method) {
  out <- args
  drop <- c("y_t", "L.x_F", "L.x_P", "L.z_F", "L.z_P", "sy", "ergm_sim_iters",
            "mu.prior_F", "Sigma.prior_F", "mu.prior_P", "Sigma.prior_P",
            "init_s_F", "init_s_P", "adapt_start", "adapt_interval",
            "gamma_F", "gamma_P", "sigma.prop_F", "sigma.prop_P")
  if (is.null(out$data) && !is.null(args$y_t)) out$data <- args$y_t
  if (is.null(out$aux.iters) && !is.null(args$ergm_sim_iters)) out$aux.iters <- args$ergm_sim_iters
  pm <- list(F = args$mu.prior_F, P = args$mu.prior_P)
  ps <- list(F = args$Sigma.prior_F, P = args$Sigma.prior_P)
  if (any(!vapply(pm, is.null, TRUE)) && is.null(out$prior.mean)) out$prior.mean <- lapply(pm, function(x) NVL(x, 0))
  if (any(!vapply(ps, is.null, TRUE)) && is.null(out$prior.sigma)) out$prior.sigma <- lapply(ps, function(x) NVL(x, 10))
  ctl <- NVL(out$control, control.bstergm())
  if (!is.null(args$adapt_start)) ctl$adapt.start <- args$adapt_start
  if (!is.null(args$adapt_interval)) ctl$adapt.interval <- args$adapt_interval
  if (!is.null(args$init_s_F) || !is.null(args$init_s_P))
    ctl$proposal.sd <- list(F = NVL(args$init_s_F, 0.025), P = NVL(args$init_s_P, 0.025))
  if (!is.null(args$gamma_F) || !is.null(args$gamma_P))
    ctl$ads.gamma <- list(F = NVL(args$gamma_F, 0.5), P = NVL(args$gamma_P, 0.5))
  if (!is.null(args$sigma.prop_F) || !is.null(args$sigma.prop_P))
    ctl$ads.var <- list(F = NVL(args$sigma.prop_F, 0.025), P = NVL(args$sigma.prop_P, 0.025))
  out$control <- ctl
  out[setdiff(names(out), drop)]
}

#' Approximate exchange algorithm with adaptive Metropolis-Hastings
#'
#' Shortcuts for [bstergm_sign()] with `method = "amh"` and `method = "ads"`.
#' They also accept the arguments of version 0.1 of the package (then called
#' \pkg{B2Lstergm}), i.e. `y_t`, `ergm_sim_iters`, `mu.prior_F`,
#' `Sigma.prior_F`, `mu.prior_P`, `Sigma.prior_P`, `init_s_F`, `init_s_P`,
#' `adapt_start`, `adapt_interval`, `gamma_F`, `gamma_P`, `sigma.prop_F` and
#' `sigma.prop_P`; the networks `L.x_F`, `L.x_P`, `L.z_F`, `L.z_P` and the
#' statistics `sy` are no longer needed, as they are computed from `y_t`.
#'
#' @param formula,data see [bstergm_sign()].
#' @param ... further arguments passed to [bstergm_sign()].
#' @return A [bstergm_sign()] fit. As in version 0.1, its elements `theta_F` and
#'   `theta_P` are arrays of posterior draws (iterations x parameters x
#'   chains).
#' @seealso [bstergm_sign()]
#' @examples
#' data(senators70_am)
#' senate <- signed_series(senators70_am[1:3])
#' fit <- aea_amh(~ edges + gwesf_pos(0.6) + gwese_pos(0.6), senate,
#'                iters = 100, burn.in = 50, aux.iters = 1000, seed = 1, verbose = FALSE)
#' dim(fit$theta_F)
#' @export
aea_amh <- function(formula, data = NULL, ...) {
  args <- .legacy_args(c(list(formula = formula, data = data), list(...)), "amh")
  fit <- do.call(bstergm_sign, c(args, list(method = "amh")))
  fit$call <- match.call()
  fit
}

#' @rdname aea_amh
#' @export
aea_ads <- function(formula, data = NULL, ...) {
  args <- .legacy_args(c(list(formula = formula, data = data), list(...)), "ads")
  fit <- do.call(bstergm_sign, c(args, list(method = "ads")))
  fit$call <- match.call()
  fit
}
