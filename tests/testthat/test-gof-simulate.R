test_that("bgof, simulate and simulate_stats work", {
  s <- small_senate()
  fit <- bstergm_sign(~ edges + gwesf_pos(0.6) + gwese_pos(0.6), s, iters = 40, burn.in = 20, seed = 1,
                 verbose = FALSE)
  gf <- bgof(fit, nsim = 5, monitor = ~ degree(0:3) + triangle_signed(), seed = 2, verbose = FALSE)
  expect_s3_class(gf, "bgof_bstergm")
  expect_equal(dim(gf$sim$F), c(5, 3 + 4 + 4, 2))
  expect_equal(dim(gf$diff$P), c(5, 11))
  expect_named(gf$groups$F, c("degree", "triangle_signed"))
  expect_output(print(gf), "p-value")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(gf))
  expect_silent(plot(gf, type = "monitor"))

  sims <- simulate(fit, nsim = 2, seed = 3)
  expect_length(sims, 2)
  expect_named(sims[[1]], c("F", "P", "y"))
  for (t in 1:2) {
    yt <- sims[[1]]$y[[t]]
    expect_true(all((yt != 0) == (s$y[[t + 1]] != 0)))
    expect_true(all((sims[[1]]$F[[t]] != 0) == (transition_networks(s)$F[[t]] != 0)))
  }
  st <- simulate(fit, nsim = 3, output = "stats", seed = 3)
  expect_equal(dim(st$F), c(3, 3, 2))

  ss <- simulate_stats(~ edges + gwesf_pos(0.6), s, coef_F = c(-1, 0.3), seed = 1)
  expect_equal(dim(ss$sz._F), c(2, 2))
  expect_null(ss$sz._P)
})
