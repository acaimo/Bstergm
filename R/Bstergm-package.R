#' Bstergm: Bayesian separable temporal ERGMs for binary and signed networks
#'
#' Bayesian inference for separable temporal exponential random graph models
#' (STERGMs). Transitions \eqn{y_{t-1} \to y_t} are split into a formation
#' network \eqn{y_{t-1} \cup y_t} and a persistence network
#' \eqn{y_{t-1} \cap y_t}, each modelled by an ERGM. For dynamic signed
#' networks (Caimo and Gollini), a signed network \eqn{y_t} is decomposed into
#' an interaction layer \eqn{x_t} (who interacts with whom) and a sign layer
#' \eqn{z_t} (which interactions are positive), and the posterior factorises
#' over the formation and persistence components of the interaction and sign
#' processes.
#'
#' @section Main functions:
#' \describe{
#'   \item{[network_series()], [signed_series()]}{prepare a sequence of binary
#'     or signed networks;}
#'   \item{[bstergm()]}{samples the posterior of a STERGM of binary networks
#'     (or of the interactions of signed networks);}
#'   \item{[bstergm_sign()]}{samples the posterior of the sign process of
#'     signed networks; [aea_amh()] and [aea_ads()] are shortcuts;}
#'   \item{[summary.bstergm()], [plot.bstergm()], [as.mcmc.list.bstergm()]}{posterior
#'     summaries, trace and density plots, and conversion to \pkg{coda};}
#'   \item{[bgof()]}{posterior predictive goodness of fit;}
#'   \item{[aux_diagnostics()]}{checks of the length of the auxiliary network
#'     simulations;}
#'   \item{[simulate.bstergm()]}{posterior predictive simulation of the observed
#'     transitions;}
#'   \item{[predict.bstergm()]}{forecasts of future binary or signed networks.}
#' }
#'
#' @section Signed-network terms:
#' The package provides native (C) ergm terms for signed networks, evaluated on
#' networks built by [as_signed_network()]:
#' [gwesf_pos][gwesf_pos-ergmTerm], [gwesf_neg][gwesf_neg-ergmTerm],
#' [gwese_pos][gwese_pos-ergmTerm], [gwese_neg][gwese_neg-ergmTerm],
#' [gwesp_signed][gwesp_signed-ergmTerm], [esp_signed][esp_signed-ergmTerm],
#' [gwdegree_neg][gwdegree_neg-ergmTerm], [degree_neg][degree_neg-ergmTerm] and
#' [triangle_signed][triangle_signed-ergmTerm], together with the sample space
#' constraint [active_signs][active_signs-ergmConstraint]. Any other ergm term
#' (e.g. `edges`, `nodematch()`, `gwdegree()`, `gwesp()`) applies to the
#' positive ties. Models written with the multilayer terms of
#' \pkg{ergm.multi} (`L()`, `gwespL()`, ...) are also supported.
#'
#' @references Caimo, A. and Gollini, I. (2025). Separable models for dynamic
#'   signed networks. arXiv:2505.07669. \doi{10.48550/arXiv.2505.07669}
#'
#' @keywords internal
"_PACKAGE"

#' @useDynLib Bstergm, .registration = TRUE
#' @import ergm
#' @import network
#' @importFrom ergm.multi Layer uncombine_network
#' @importFrom statnet.common NVL list_rhs.formula
#' @importFrom stats coef simulate predict rnorm runif sd quantile var cov median
#'   setNames glm.fit binomial vcov
#' @importFrom graphics par plot lines abline hist boxplot points legend
#'   matplot axis mtext title
#' @importFrom grDevices adjustcolor devAskNewPage
#' @importFrom utils head tail
#' @importFrom methods is
NULL
