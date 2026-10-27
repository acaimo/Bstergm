test_that("init = 'ads' starts the AMH chains from an initial ADS population", {
  s <- small_senate()
  f <- ~ edges + gwesf_pos(0.6) + gwese_neg(0.6)
  ctl <- control.bstergm(ads.init.iters = 20)
  fit <- bstergm_sign(f, s, iters = 30, burn.in = 10, n_chains = 2, init = "ads", seed = 1,
                 verbose = FALSE, control = ctl)
  expect_equal(dim(fit$theta$F), c(20, 3, 2))
  expect_named(fit$init.ads, c("F", "P"))
  expect_equal(fit$init.ads$F$chains, 6L)
  expect_equal(dim(fit$init.ads$F$cov), c(3, 3))
  expect_true(all(is.finite(fit$proposal$F[[1]])))
  if (.Platform$OS.type == "unix") {
    fit2 <- bstergm_sign(f, s, iters = 30, burn.in = 10, n_chains = 2, init = "ads", seed = 1,
                    verbose = FALSE, control = ctl, cores = 2)
    expect_identical(fit$theta, fit2$theta)
  }
  expect_error(bstergm_sign(f, s, method = "ads", init = "ads", iters = 5, burn.in = 1, verbose = FALSE),
               "only available")
})

test_that("predict() forecasts valid signed networks, independently of the cores", {
  s <- small_senate(20)
  fz <- bstergm_sign(~ edges + gwesf_pos(0.6), s, iters = 20, burn.in = 10, seed = 1, verbose = FALSE,
                aux.iters = 300)
  fx <- bstergm(~ edges + gwesp(0.6, fixed = TRUE), s, iters = 20,
                burn.in = 10, seed = 1, verbose = FALSE, aux.iters = 300)
  fc <- predict(fz, fx, horizon = 2, nsim = 3, visits = 2, seed = 4, verbose = FALSE)
  expect_s3_class(fc, "bstergm_forecast")
  expect_length(fc$sims, 3)
  expect_named(fc$sims[[1]], c("T+1", "T+2"))
  for (tr in fc$sims) for (y in tr) {
    expect_true(isSymmetric(unname(y)))
    expect_true(all(y %in% c(-1L, 0L, 1L)))
  }
  sm <- summary(fc)
  expect_equal(dim(sm[["T+1"]]), c(7, 4))
  expect_output(print(fc), "Forecast")
  if (.Platform$OS.type == "unix") {
    fc2 <- predict(fx, fz, horizon = 2, nsim = 3, visits = 2, seed = 4, verbose = FALSE, cores = 2)
    expect_identical(fc$sims, fc2$sims)
  }
  expect_error(predict(fz, fz, verbose = FALSE), "one bstergm_sign\\(\\) fit")
})

test_that("posterior predictive simulations do not depend on the number of cores", {
  skip_on_os("windows")
  s <- small_senate(20)
  fit <- bstergm_sign(~ edges + gwese_pos(0.6), s, iters = 20, burn.in = 10, seed = 1, verbose = FALSE,
                 aux.iters = 300)
  g1 <- bgof(fit, nsim = 4, seed = 2, verbose = FALSE)
  g2 <- bgof(fit, nsim = 4, seed = 2, verbose = FALSE, cores = 2)
  expect_identical(g1$sim, g2$sim)
  s1 <- simulate(fit, nsim = 3, seed = 5)
  s2 <- simulate(fit, nsim = 3, seed = 5, cores = 2)
  expect_identical(s1, s2)
})
