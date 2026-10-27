## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4.5,
                      fig.align = "center", out.width = "100%")

## ----data---------------------------------------------------------------------
library(Bstergm)
data(senators70_am)

senate <- signed_series(senators70_am[1:3],
                        vertex.attr = attr(senators70_am$S99, "senators"))
senate

## ----summary------------------------------------------------------------------
summary(senate)

## ----plot-data, fig.height = 3------------------------------------------------
plot(senate, vertex.col = "Party")

## ----stats--------------------------------------------------------------------
nw <- as_signed_network(senators70_am$S99,
                        vertex.attr = attr(senators70_am$S99, "senators"))
summary(nw ~ edges + nodematch("Party", diff = TRUE, levels = "R") +
          gwesf_pos(0.6) + gwese_pos(0.6) + gwese_neg(0.6) + triangle_signed())

## ----multilayer---------------------------------------------------------------
positive <- ~positive
negative <- ~(active & !positive)
z <- network::as.network((senators70_am$S99 > 0) * 1, directed = FALSE)
x <- network::as.network((senators70_am$S99 != 0) * 1, directed = FALSE)
summary(~ gwespL(0.6, fixed = TRUE, L.base = positive, Ls.path = list(positive, positive)) +
          gwespL(0.6, fixed = TRUE, L.base = positive, Ls.path = list(negative, negative)) +
          gwespL(0.6, fixed = TRUE, L.base = negative, Ls.path = list(negative, negative)),
        basis = Layer(positive = z, active = x))

## ----quick--------------------------------------------------------------------
f <- ~ edges + nodematch("Party", diff = TRUE, levels = "R") +
  gwdegree(0.2, fixed = TRUE) + gwesf_pos(0.6) + gwese_pos(0.6) + gwese_neg(0.6)

quick <- bstergm_sign(f, data = senate, iters = 300, burn.in = 100, seed = 1, verbose = FALSE)
quick

## ----load, include = FALSE----------------------------------------------------
fit <- readRDS(system.file("extdata", "senate_fit.rds", package = "Bstergm"))
fit_x <- readRDS(system.file("extdata", "senate_fit_interaction.rds", package = "Bstergm"))
mins <- round(fit$time$total / 60)

## ----long, eval = FALSE-------------------------------------------------------
# fit <- bstergm_sign(f, data = senate, iters = 35000, burn.in = 10000, thin = 25,
#                     n_chains = 2, cores = 2, seed = 2025,
#                     prior.mean = c(-1, 0, 0, 0, 0, 0), prior.sigma = 25)

## ----load-show, eval = FALSE--------------------------------------------------
# fit <- readRDS(system.file("extdata", "senate_fit.rds", package = "Bstergm"))

## ----summary-fit--------------------------------------------------------------
summary(fit)

## ----trace, fig.height = 9----------------------------------------------------
plot(fit, components = "F")

## ----interpret----------------------------------------------------------------
sapply(fit$theta, function(a) round(apply(a, 2, function(v) mean(v > 0)), 3))

## ----gof----------------------------------------------------------------------
gf <- bgof(fit, nsim = 50, seed = 1, verbose = FALSE,
           monitor = ~ degree(0:12) + degree_neg(0:12) + triangle_signed())
gs <- summary(gf)
gs$F[c(1:6, 33:36), ]
plot(gf, type = "model")

## ----gof-monitor, fig.height = 7----------------------------------------------
plot(gf, type = "monitor")

## ----auxdiag, fig.height = 5.5------------------------------------------------
dg <- aux_diagnostics(fit, visits = 50, seed = 1)
dg
plot(dg)

## ----auxdiag-values, include = FALSE------------------------------------------
trF <- dg$results$F$traj$t1[[1]]$stats
half <- seq_len(nrow(trF)) > nrow(trF) / 2
lr <- round(colMeans(trF[half, ]))

## ----interaction, eval = FALSE------------------------------------------------
# fit_x <- bstergm(~ edges + nodematch("Party", diff = TRUE, levels = "R") +
#                    gwdegree(0.2, fixed = TRUE) + gwesp(0.6, fixed = TRUE),
#                  data = senate, iters = 35000,
#                  burn.in = 10000, thin = 25, n_chains = 2, cores = 2, seed = 2025,
#                  prior.mean = c(-1, 0, 0, 0), prior.sigma = 25)

## ----interaction-summary------------------------------------------------------
summary(fit_x)

## ----interaction-ridge--------------------------------------------------------
round(cor(as.matrix(as.mcmc.list(fit_x, components = "F"))), 2)

## ----simulate-----------------------------------------------------------------
sims <- simulate(fit, nsim = 1, seed = 2)
sapply(sims[[1]]$y, signed_triad_census)

## ----predict-binary-----------------------------------------------------------
summary(predict(fit_x, horizon = 1, nsim = 20, seed = 3, verbose = FALSE))

## ----predict------------------------------------------------------------------
fc <- predict(fit, fit_x, horizon = 1, nsim = 20, seed = 3, verbose = FALSE)
summary(fc)

