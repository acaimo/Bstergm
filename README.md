
# **Bstergm**: Bayesian inference for separable temporal
exponential random graph models

<!-- badges: start -->
<!-- badges: end -->

**Bstergm** implements Bayesian inference for separable temporal
exponential random graph models (STERGMs) of dynamic **binary networks**,
with `bstergm()`, and of dynamic **signed networks**, with `bstergm_sign()`,
following Caimo and Gollini (2025). It is the successor of **B2Lstergm** and
**Bstergm.sign**.

Each transition $y_{t-1} \to y_t$ is split into the *formation* network
$y_{t-1} \cup y_t$ and the *persistence* network $y_{t-1} \cap y_t$, each
modelled by an ERGM, so the posterior factorises into independent formation
and persistence components. A signed network is further decomposed into an
*interaction* layer $x_t$ (who interacts with whom), a binary network fitted
with `bstergm()`, and a *sign* layer $z_t$ (which interactions are positive),
fitted with `bstergm_sign()`. The doubly intractable posteriors are sampled
with the approximate exchange algorithm, using adaptive Metropolis-Hastings
or adaptive direction sampling.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("acaimo/Bstergm")
```

The package contains C code linked against **ergm**, so a compiler is needed
when installing from source (Rtools on Windows).

## Example

``` r
library(Bstergm)
data(senators70_am)

# Signed networks of 70 U.S. senators in the 99th, 100th and 101st Congress
senate <- signed_series(senators70_am[1:3],
                        vertex.attr = attr(senators70_am$S99, "senators"))
summary(senate)

# Sign process: positive ties, Republican homophily, degree, and balance
# (edgewise shared friends/enemies)
fit <- bstergm_sign(~ edges + nodematch("Party", diff = TRUE, levels = "R") +
                      gwdegree(0.2, fixed = TRUE) +
                      gwesf_pos(0.6) + gwese_pos(0.6) + gwese_neg(0.6),
                    data = senate, iters = 35000, burn.in = 10000, thin = 25,
                    n_chains = 2, cores = 4, seed = 1)
summary(fit)
plot(fit)

# Posterior predictive goodness of fit, with extra monitored statistics,
# and checks of the auxiliary simulations
gf <- bgof(fit, nsim = 100, monitor = ~ degree(0:10) + triangle_signed())
plot(gf)
dg <- aux_diagnostics(fit)
plot(dg)

# Binary networks: the interactions, with ordinary ergm terms. Any list of
# binary adjacency matrices or networks can be used (see network_series())
fit_x <- bstergm(~ edges + nodematch("Party", diff = TRUE, levels = "R"),
                 data = senate, cores = 4)

# Forecast the next binary network, or the next signed network from both fits
summary(predict(fit_x, horizon = 1, nsim = 100))
summary(predict(fit, fit_x, horizon = 1, nsim = 100))
```

All plots are in grayscale. `method = "ads"` (or `aea_ads()`) samples with
adaptive direction sampling, and `init = "ads"` runs adaptive direction sampling to find starting values
before adaptive Metropolis-Hastings, as in the paper.

## Signed-network terms

The model formula of `bstergm_sign()` is evaluated on the positive ties of
each formation or persistence network, with the interactions held fixed.
Ordinary ergm terms (`edges`, `nodematch()`, `gwdegree()`, `edgecov()`, ...)
therefore refer to positive ties. The terms below are implemented natively in
C:

| term | statistic |
|---|---|
| `gwesf_pos(decay)` | positive ties weighted by their shared friends ($\mathrm{gwesf}^+$, balanced $+\,+\,+$) |
| `gwese_pos(decay)` | positive ties weighted by their shared enemies ($\mathrm{gwese}^+$, balanced $+\,-\,-$) |
| `gwesf_neg(decay)` | negative ties weighted by their shared friends ($\mathrm{gwesf}^-$, unbalanced $+\,+\,-$) |
| `gwese_neg(decay)` | negative ties weighted by their shared enemies ($\mathrm{gwese}^-$, $-\,-\,-$) |
| `gwesp_signed(decay, base, path)` | general version, also with mixed $+/-$ paths |
| `esp_signed(d, base, path)` | edgewise shared partner distribution by sign |
| `gwdegree_neg(decay)`, `degree_neg(d)` | degree of negative ties |
| `triangle_signed(type)` | signed triangle census |

Multilayer formulas of **ergm.multi** (`L()`, `gwespL()`, ...) on the layers
`positive` and `active`, as used in version 0.1, are also accepted. When
possible they are translated automatically into the equivalent native terms.

## Speed

Seconds per iteration of the approximate exchange algorithm for the model of
the paper on the U.S. Senate data (6 parameters, formation and persistence
components, 2 transitions), measured on the same two-core machine:

| implementation | auxiliary steps | s / iteration | speed-up |
|---|---|---:|---:|
| version 0.1 (B2Lstergm) | 5000 | 3.67 | 1 |
| this version, `ergm.multi` formula | 5000 | 0.318 | 12 |
| this version, native signed terms | 5000 | 0.051 | 72 |
| this version, native terms, default length | 3 per free dyad | 0.037 | 99 |
| this version, native terms, default length, 2 cores | 3 per free dyad | 0.032 | 113 |

The analysis of the vignette (two chains of 35,000 iterations) takes 19
minutes on two cores, against about 3 days with version 0.1. The gain comes
from:

- building the ergm model, proposal and constraints of each transition once,
  instead of at every iteration (version 0.1 spent about 3.3 s per iteration
  on this set-up, against about 0.4 s of actual simulation);
- native C signed-network statistics with shared-partner caches: 2.2
  microseconds per proposal and 0.8 ms of set-up per simulation, against 19.6
  microseconds and 12.4 ms with the multilayer terms of **ergm.multi**;
- proposing only sign changes of interacting dyads: in version 0.1 half of the
  proposals toggled the fixed interaction layer and were always rejected;
- running the formation and persistence components, and the chains, in
  parallel.

## Checking the approximation

The approximate exchange algorithm treats the end of each auxiliary
simulation, started at the observed network, as a draw from the model. Use
`aux_diagnostics()` and `bgof()` to check this. For near-degenerate models,
long simulations drift far from the observed networks and the sampled
distribution depends on the length of the auxiliary simulations
(`aux.iters`, `control.bstergm(aux.visits =)`). The vignette works through
this for the model of the paper.

## Coming from Bstergm.sign or B2Lstergm 0.1

Replace `library(Bstergm.sign)` or `library(B2Lstergm)` by `library(Bstergm)`.
Fits of signed networks, `bstergm(...)` in Bstergm.sign, are now
`bstergm_sign(...)`, and `bstergm(..., process = "interaction")` is now
`bstergm(...)`. The option `Bstergm.sign.cores` is now `Bstergm.cores`. `aea_amh()`,
`aea_ads()` and `simulate_stats()` accept the arguments of version 0.1, and
the networks `L.x_F`, `L.x_P`, `L.z_F`, `L.z_P` and statistics `sy` are now
computed internally from `y_t`. Version 0.1 had bugs that affected its
results; see `NEWS.md`.

## Reference

Caimo, A., & Gollini, I. (2026). Separable models for dynamic signed networks. 
Journal of the Royal Statistical Society Series A: Statistics in Society, qnag053.
