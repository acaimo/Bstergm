test_that("bstergm returns well-formed, reproducible output", {
  s <- small_senate()
  f <- ~ edges + gwesf_pos(0.6) + gwese_pos(0.6)
  fit <- bstergm_sign(f, s, iters = 60, burn.in = 20, n_chains = 2, seed = 1, verbose = FALSE)
  expect_s3_class(fit, "bstergm")
  expect_equal(dim(fit$theta$F), c(40, 3, 2))
  expect_identical(fit$theta_F, fit$theta$F)
  expect_equal(fit$engine, "native")
  expect_true(all(fit$acceptance$F >= 0 & fit$acceptance$F <= 1))
  fit2 <- bstergm_sign(f, s, iters = 60, burn.in = 20, n_chains = 2, seed = 1, verbose = FALSE)
  expect_identical(fit$theta, fit2$theta)
  if (.Platform$OS.type == "unix") {
    fit3 <- bstergm_sign(f, s, iters = 60, burn.in = 20, n_chains = 2, seed = 1, cores = 2, verbose = FALSE)
    expect_identical(fit$theta, fit3$theta)
  }
  sm <- summary(fit)
  expect_s3_class(sm, "summary.bstergm")
  expect_true(all(c("Mean", "SD", "ESS", "Rhat") %in% colnames(sm$tables$F)))
  expect_output(print(fit), "Formation")
  expect_length(coef(fit)$P, 3)
  ml <- as.mcmc.list(fit)
  expect_equal(coda::nchain(ml), 2)
  expect_equal(coda::nvar(ml), 6)
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(fit, ask = FALSE))
})

test_that("thinning, ADS and the binary STERGM work", {
  s <- small_senate()
  fit <- bstergm_sign(~ edges + gwese_neg(0.6), s, components = "P", method = "ads", n_chains = 4,
                 iters = 30, burn.in = 10, thin = 5, seed = 2, verbose = FALSE)
  expect_equal(dim(fit$theta$P), c(4, 2, 4))
  expect_null(fit$theta$F)
  fi <- bstergm(~ edges + gwesp(0.5, fixed = TRUE), s, iters = 30, burn.in = 10, seed = 3, verbose = FALSE)
  expect_equal(fi$engine, "standard")
  expect_error(bstergm(~ edges + gwese_pos(0.5), s, verbose = FALSE),
               "need bstergm_sign")
})

test_that("the multilayer formulation of version 0.1 is supported", {
  s <- small_senate(20)
  positive <- ~positive
  negative <- ~(active & !positive)
  f <- ~ L(~edges, ~active) + L(~edges, positive) +
    gwespL(0.6, fixed = TRUE, L.base = positive, Ls.path = list(negative, negative)) +
    gwespL(0.6, fixed = TRUE, L.base = negative, Ls.path = list(positive, negative))
  fit <- bstergm_sign(f, s, iters = 15, burn.in = 5, seed = 4, verbose = FALSE, aux.iters = 500)
  expect_equal(fit$engine, "native")
  expect_equal(dimnames(fit$theta$F)[[2]],
               c("edges", "gwese_pos.fixed.0.6", "gwesp_signed(-,+-).fixed.0.6"))
  expect_identical(fit$formula.original, f)
  fl <- suppressWarnings(bstergm_sign(f, s, iters = 15, burn.in = 5, seed = 4, verbose = FALSE,
                                 aux.iters = 500, control = control.bstergm(translate = FALSE)))
  expect_equal(fl$engine, "layer")
  expect_named(fl$fixed$F, "L(active)~edges")
  expect_equal(dim(fl$theta$F)[2], 3)
  expect_equal(unname(fl$obs.stats$F[, -1]), unname(fit$obs.stats$F), tolerance = 1e-10)
  f2 <- ~ L(~edges, positive) + L(~triangles, ~active & positive)
  expect_null(Bstergm:::.translate_formula(f2))
})

test_that("the exchange algorithm recovers an exact (dyad-independent) posterior", {
  skip_on_cran()
  s <- small_senate(40)
  fit <- bstergm_sign(~ edges + nodematch("Party"), s, components = "F", iters = 2500, burn.in = 500,
                 seed = 5, verbose = FALSE, prior.sigma = 10, control = control.bstergm(aux.visits = 5))
  tn <- transition_networks(s)
  same <- outer(s$vertex.attr$Party, s$vertex.attr$Party, "==")
  X <- NULL; y <- NULL
  for (t in 1:2) {
    m <- tn$F[[t]]; free <- upper.tri(m) & m != 0
    X <- rbind(X, cbind(1, same[free])); y <- c(y, m[free] > 0)
  }
  b <- c(0, 0)
  for (i in 1:50) {
    mu <- plogis(drop(X %*% b))
    H <- crossprod(X, X * (mu * (1 - mu))) + diag(0.1, 2)
    b <- b + solve(H, crossprod(X, y - mu) - 0.1 * b)
  }
  sm <- summary(fit)$tables$F
  expect_lt(max(abs(unname(sm[, "Mean"]) - drop(b))), 0.06)
  expect_lt(max(abs(unname(sm[, "SD"]) / sqrt(diag(solve(H))) - 1)), 0.3)
})

test_that("arguments are validated", {
  s <- small_senate(10)
  expect_error(bstergm_sign(~ edges, s, iters = 10, burn.in = 10, verbose = FALSE), "smaller")
  expect_error(bstergm_sign(~ edges, s, method = "ads", n_chains = 2, iters = 5, burn.in = 1, verbose = FALSE),
               "at least 3|>= 3")
  expect_error(bstergm_sign(~ edges, s, prior.mean = 1:3, iters = 5, burn.in = 1, verbose = FALSE), "length")
  expect_error(bstergm_sign(list(~edges, ~edges, ~edges), s, iters = 5, burn.in = 1, verbose = FALSE),
               "one formula per transition")
})

test_that("a fit can be continued and diagnosed", {
  s <- small_senate()
  f <- ~ edges + gwesf_pos(0.6)
  fit <- bstergm_sign(f, s, iters = 30, burn.in = 10, n_chains = 2, seed = 1, verbose = FALSE)
  fit2 <- bstergm_sign(f, s, iters = 20, burn.in = 0, n_chains = 2, init = fit, seed = 2, verbose = FALSE)
  expect_equal(dim(fit2$theta$F), c(20, 2, 2))
  expect_error(bstergm_sign(~ edges, s, init = fit, iters = 5, burn.in = 1, verbose = FALSE), "does not match")
  dg <- aux_diagnostics(fit, visits = 4, nsamples = 8, seed = 1)
  expect_s3_class(dg, "aux_diagnostics_bstergm")
  expect_equal(dim(dg$results$F$traj$t1[[1]]$stats), c(9, 2))
  expect_output(print(dg), "drift")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(dg))
})
