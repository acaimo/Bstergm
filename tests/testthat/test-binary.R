test_that("bstergm fits binary STERGMs from matrices, networks and signed data", {
  s <- small_senate(25)
  bin <- lapply(s$y, function(m) abs(m))
  ns <- network_series(bin, vertex.attr = s$vertex.attr)
  expect_s3_class(ns, "network_series")
  expect_output(print(ns), "binary networks")
  f <- ~ edges + nodematch("Party")
  fit <- bstergm(f, ns, iters = 30, burn.in = 10, aux.iters = 300, seed = 1, verbose = FALSE)
  expect_s3_class(fit, "bstergm")
  expect_equal(fit$process, "interaction")
  expect_equal(dim(fit$theta$F), c(20, 2, 1))
  fit2 <- bstergm(f, s, iters = 30, burn.in = 10, aux.iters = 300, seed = 1, verbose = FALSE)
  expect_identical(fit$theta, fit2$theta)
  nets <- lapply(bin, function(m) {
    nw <- network::as.network(m, directed = FALSE)
    network::set.vertex.attribute(nw, "Party", s$vertex.attr$Party)
  })
  fit3 <- bstergm(f, nets, iters = 30, burn.in = 10, aux.iters = 300, seed = 1, verbose = FALSE)
  expect_identical(unname(fit$theta$F), unname(fit3$theta$F))
  expect_output(print(fit), "binary network")
  fc <- predict(fit, horizon = 2, nsim = 2, visits = 2, seed = 1, verbose = FALSE)
  expect_false(fc$signed)
  expect_true(all(unlist(fc$sims) %in% c(0L, 1L)))
  expect_equal(rownames(summary(fc)[["T+1"]]), c("ties", "triangles"))
  expect_error(predict(bstergm_sign(~ edges, s, iters = 5, burn.in = 1, aux.iters = 100, verbose = FALSE)),
               "also needs a bstergm")
  expect_error(bstergm(f, list(matrix(0, 3, 3))), "two time points")
})

test_that("fits do not reference the package namespace", {
  ctrl <- control.bstergm()
  expect_identical(environment(ctrl$MCMC.prop), globalenv())
  for (f in c("senate_fit.rds", "senate_fit_interaction.rds")) {
    x <- readRDS(system.file("extdata", f, package = "Bstergm"))
    expect_false(grepl("Bstergm", rawToChar(serialize(x, NULL, ascii = TRUE)), fixed = TRUE))
  }
})
