# Bstergm 0.2.0

The package was renamed from **B2Lstergm** (and from the intermediate name
**Bstergm.sign**) to **Bstergm**. It now covers dynamic binary networks with
`bstergm()` and dynamic signed networks with `bstergm_sign()`. Its
computational core was rewritten, several bugs that affected the results of
version 0.1 were fixed, and the package was completed with data handling,
native signed-network statistics, posterior summaries, goodness of fit,
diagnostics, simulation and forecasting. All plots are in grayscale.

Scripts written for version 0.1 keep working after replacing
`library(B2Lstergm)` by `library(Bstergm)`: `aea_amh()`, `aea_ads()` and
`simulate_stats()` accept their old arguments, `Layer()` is re-exported, and
multilayer formulas with the constraint term `L(~edges, active)` are handled
automatically.

Coming from **Bstergm.sign**: `bstergm()` for signed networks is now
`bstergm_sign()`; `bstergm(..., process = "interaction")` is now
`bstergm()`; the option `Bstergm.sign.cores` is now `Bstergm.cores`.

## Bug fixes (these change the results of version 0.1)

* `simulate_stats()`, used at every iteration of `aea_amh()` and `aea_ads()`:
  from the second transition onwards, each simulation started from the
  network simulated for the previous transition. Its positive ties outside
  the current interactions were then kept fixed, giving positive ties without
  an interaction and wrong statistics (e.g. 145 of 856 positive ties in one
  simulation of the second transition of the vignette data, at zero
  parameters). With two or more transitions, as in the vignette and the
  paper, the posterior was therefore biased. Each transition now starts from
  its observed network.
* `aea_ads()`: the adaptive direction sampling move used the difference of two
  scalar entries of the flattened parameter matrix, added to every parameter,
  instead of the difference between the parameter vectors of two other chains.
  These entries could be parameters of the chain being updated (e.g. for the
  first two chains with 10 chains and 6 parameters), in which case the move
  was not symmetric and the sampler did not target the posterior.
* The reconstruction of the network at time t used
  `y_P | (y_F - y_{t-1})`, in which negative differences count as `TRUE`.
  It now uses `y_P | (y_F & !y_{t-1})`, restricted to the observed
  interactions.
* Vignette: the covariate `rep` was 1 in every row of a Republican senator.
  `edgecov(rep)` therefore counted positive ties whose lower-indexed senator is
  a Republican, not ties between two Republicans. The vignette now uses
  `nodematch("Party", diff = TRUE, levels = "R")`.
* `aea_amh()` and `aea_ads()` failed with fewer than 100 iterations. Progress
  reporting and some storage were hard-wired to 50,000 iterations.

## Computational improvements

On the U.S. Senate data, one iteration of the model of the paper (formation
and persistence, two transitions) takes 0.037 s with the default auxiliary
length (0.051 s with 5,000 auxiliary steps), against 3.67 s for version 0.1
with 5,000 auxiliary steps: about 100 times faster on one core. The analysis
of the vignette (two chains of 35,000 iterations) takes 19 minutes on two
cores, against about 3 days with version 0.1.

* The ergm model, MCMC proposal and sample-space constraints of every
  transition are built once, before sampling. Version 0.1 rebuilt them with
  `simulate()` and `summary()` at every iteration: about 3.3 s per iteration,
  against about 0.4 s of actual simulation.
* Native C implementations of the signed-network statistics
  (`gwesf_pos()`, `gwese_pos()`, `gwese_neg()`, ...). They keep caches of
  shared partners by sign and cost O(degree) per proposal: 2.2 microseconds
  per proposal on the Senate data, against 19.6 microseconds with the
  equivalent multilayer terms of **ergm.multi**, which give the same values
  and, with the same seed, the same simulations.
* The auxiliary simulations propose only sign changes of interacting dyads.
  In version 0.1 the dyads of the fixed interaction layer were also proposed,
  and always rejected through the `L(~edges, active)` term with an infinite
  coefficient: half of the proposals were wasted.
* The sampler takes the change in the statistics directly from the C
  sampler; no `summary()` calls are needed.
* The formation and persistence components (and, for adaptive
  Metropolis-Hastings, the chains) are independent and run in parallel with
  `cores =`. Every chain and every simulation has its own L'Ecuyer-CMRG
  random number stream, so results are identical for any number of cores.
* Chains start at a penalised maximum pseudo-likelihood estimate pooled over
  transitions, whose covariance also scales the initial proposal. This
  shortens burn-in considerably compared with starting at the prior mean.
* The length of the auxiliary simulations scales with the number of free
  dyads of each transition (`control.bstergm(aux.visits =)`).
* Adaptive Metropolis-Hastings updates the proposal covariance recursively
  instead of recomputing it from the stored history.

## New features

* `signed_series()` builds the data from a list of signed adjacency matrices,
  signed `network` objects or the two-layer `Layer()` networks of version
  0.1. `transition_networks()` returns the formation and persistence networks,
  so the inputs `L.x_F`, `L.x_P`, `L.z_F`, `L.z_P` and `sy` no longer need
  to be prepared by hand. `summary()` and `plot()` methods describe the data.
* `bstergm()` fits separable temporal ERGMs of binary networks with ordinary
  ergm terms; `bstergm_sign()` fits the sign process of signed networks.
  `network_series()` builds a binary series from adjacency matrices or
  `network` objects; a signed series can also be passed to `bstergm()`, which
  then models its interactions. Both fit the formation and persistence
  components with `method = "amh"` or `method = "ads"`, with
  component-specific priors, starting values and offsets. `aea_amh()` and
  `aea_ads()` remain as shortcuts of `bstergm_sign()`.
* `init = "ads"` runs adaptive direction sampling first and starts the
  adaptive Metropolis-Hastings chains from its final population, the scheme
  of Caimo and Gollini. A previous fit can be passed as `init` to continue
  its chains.
* Signed-network ergm terms: `gwesf_pos`, `gwesf_neg`, `gwese_pos`,
  `gwese_neg`, `gwesp_signed`, `esp_signed`, `gwdegree_neg`, `degree_neg`
  and `triangle_signed`, plus the constraint `active_signs`. With
  `as_signed_network()` they also make it possible to fit static signed
  ERGMs with `ergm::ergm()`.
* Models written with the multilayer terms of **ergm.multi** (as in version
  0.1) are still supported. They are translated automatically into the
  equivalent native terms; the translation is used only if it reproduces the
  observed statistics exactly. Formulas that cannot be translated run with the
  multilayer terms. Statistics that are constant on the sample space, such as
  `L(~edges, ~active)`, are detected and dropped.
* Time-varying covariates, by passing a list of formulas (one per
  transition). Lagged dyadic covariates of the previous time point are
  available as `edgecov("lag.pos")`, `edgecov("lag.neg")` and
  `edgecov("lag.active")`.
* `summary()` (posterior means, standard deviations, quantiles, effective
  sample sizes and Gelman-Rubin diagnostics), `print()`, `coef()`, `plot()`
  (densities, traces, autocorrelations) and `as.mcmc.list()` / `as.mcmc()`
  for **coda**.
* `bgof()`: posterior predictive goodness of fit with the centred cumulative
  differences of the paper, plus any monitored statistics such as degree
  distributions or the signed triangle census. Its simulations are five times
  longer than the auxiliary simulations used for fitting, so that the check is
  not tied to the fitting length.
* `aux_diagnostics()` runs long simulations from the observed networks at
  posterior values, to check whether the auxiliary simulations are long
  enough and whether the model is near-degenerate there.
* `simulate()` produces posterior predictive formation, persistence and
  full networks of the observed transitions. `predict()` forecasts future
  binary networks from a `bstergm()` fit, and future signed networks from a
  `bstergm_sign()` fit together with a `bstergm()` fit of the interactions.
* `signed_triad_census()`, and a test suite, including a check that the
  sampler recovers an exactly known posterior.

## Notes on the model of the paper

* The approximate exchange algorithm treats the end of each auxiliary
  simulation, started at the observed network, as an exact draw from the
  model. For the sign model of the paper on the Senate data, long simulations
  at the posterior means move far away from the observed networks: the
  formation networks go from 26-30% positive ties to about 45% at the
  estimates of this version and about 95% at those reported in the paper. The
  sampled distribution is then a local, contrastive-divergence-like
  approximation that depends on the length of the auxiliary simulations.
  `aux_diagnostics()` and `bgof()` reveal this; the vignette shows how.
* In the model of the interactions (a `bstergm()` fit), the change statistics of `gwdegree(0.2)` and
  (for formation) `gwesp(0.6)` are zero or almost constant on these dense
  networks, so these parameters are not identified separately from `edges`;
  see the vignette.

# B2Lstergm 0.1.0

* First version: `aea_amh()`, `aea_ads()`, `simulate_stats()` and the
  `senators70_am` data.
