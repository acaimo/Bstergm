test_that("native signed statistics agree with brute force", {
  set.seed(1)
  for (rep in 1:10) {
    m <- random_signed(sample(6:25, 1), runif(1, 0.2, 0.9), runif(1, 0.2, 0.8))
    nw <- as_signed_network(m)
    expect_equal(unname(summary(native_formula, basis = nw)), unname(brute_signed_stats(m)),
                 tolerance = 1e-10)
    expect_equal(unname(summary(nw ~ triangle_signed())), unname(signed_triad_census(m)))
  }
})

test_that("change statistics are consistent along simulations, also off the active set", {
  set.seed(2)
  m <- random_signed(15, 0.5)
  nw <- as_signed_network(m)
  sim <- simulate(native_formula, basis = nw, coef = rnorm(9, 0, 0.2), nsim = 3,
                  control = control.simulate.formula(MCMC.burnin = 300, MCMC.interval = 50),
                  output = "network")
  for (s in sim) {
    expect_equal(unname(summary(native_formula, basis = s)),
                 unname(brute_signed_stats(signed_adjacency(s))), tolerance = 1e-10)
  }
  st <- attr(sim, "stats")
  expect_equal(unname(st[3, ]), unname(brute_signed_stats(signed_adjacency(sim[[3]]))), tolerance = 1e-8)
})

test_that("native terms reproduce ergm.multi's multilayer statistics", {
  set.seed(3)
  m <- random_signed(12, 0.6)
  z <- network::as.network((m > 0) * 1, directed = FALSE)
  x <- network::as.network((m != 0) * 1, directed = FALSE)
  L <- ergm.multi::Layer(positive = z, active = x)
  pos <- ~positive
  neg <- ~(active & !positive)
  fL <- ~ gwespL(0.6, fixed = TRUE, L.base = pos, Ls.path = list(pos, pos)) +
    gwespL(0.6, fixed = TRUE, L.base = pos, Ls.path = list(neg, neg)) +
    gwespL(0.6, fixed = TRUE, L.base = neg, Ls.path = list(neg, neg))
  expect_equal(unname(summary(fL, basis = L)),
               unname(summary(as_signed_network(m) ~ gwesf_pos(0.6) + gwese_pos(0.6) + gwese_neg(0.6))),
               tolerance = 1e-10)
})

test_that("active_signs keeps the interactions fixed", {
  set.seed(4)
  m <- random_signed(20, 0.5)
  nw <- as_signed_network(m)
  sim <- simulate(nw ~ edges + gwese_pos(0.6), coef = c(0.3, 0.2), constraints = ~active_signs,
                  control = control.simulate.formula(MCMC.burnin = 2000))
  expect_identical(unname(signed_adjacency(sim) != 0), m != 0)
})

test_that("signed network conversions round-trip", {
  set.seed(5)
  m <- random_signed(10)
  expect_identical(unname(signed_adjacency(as_signed_network(m))), unname(m))
  z <- network::as.network((m > 0) * 1, directed = FALSE)
  x <- network::as.network((m != 0) * 1, directed = FALSE)
  expect_identical(unname(signed_adjacency(ergm.multi::Layer(positive = z, active = x))), unname(m))
  expect_error(signed_adjacency(matrix(2, 3, 3)), "only contain")
  expect_error(as_signed_network(network::network.initialize(3, directed = FALSE)), "edge attribute")
  expect_error(summary(network::network.initialize(4, directed = FALSE) ~ gwese_pos(0.5)), "active")
})
