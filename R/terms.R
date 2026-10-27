.active_el <- function(nw) {
  el <- nw %n% ".active"
  if (is.null(el)) {
    stop("Signed-network terms need the active (interacting) dyads in the network attribute ",
         sQuote(".active"), ". Build the network with as_signed_network().", call. = FALSE)
  }
  if (is.network(el)) el <- as.matrix(el, matrix.type = "edgelist")
  el <- as.matrix(el)
  if (length(el) == 0) return(matrix(integer(0), 0, 2))
  if (ncol(el) != 2) stop("Network attribute '.active' must be a two-column edge list.", call. = FALSE)
  n <- network.size(nw)
  el <- cbind(pmin(el[, 1], el[, 2]), pmax(el[, 1], el[, 2]))
  storage.mode(el) <- "integer"
  if (any(el < 1L | el > n)) stop("Network attribute '.active' refers to non-existent vertices.", call. = FALSE)
  el <- el[el[, 1] != el[, 2], , drop = FALSE]
  unique(el)
}

.active_adj <- function(nw, el = .active_el(nw)) {
  n <- network.size(nw)
  Matrix::sparseMatrix(i = c(el[, 1], el[, 2]), j = c(el[, 2], el[, 1]),
                       x = 1, dims = c(n, n))
}

.parse_base <- function(base) {
  base <- as.character(base)[1]
  if (base %in% c("+", "1", "pos", "positive")) return(1L)
  if (base %in% c("-", "-1", "neg", "negative")) return(-1L)
  stop("Argument 'base' must be \"+\" or \"-\".", call. = FALSE)
}

.parse_path <- function(path) {
  path <- gsub("\\s", "", as.character(path)[1])
  switch(path,
         "++" = 0L, "--" = 1L, "+-" = 2L, "-+" = 2L,
         stop("Argument 'path' must be one of \"++\", \"--\" or \"+-\".", call. = FALSE))
}

.gw_table <- function(decay, len) exp(decay) * (1 - (1 - exp(-decay))^(seq_len(len) - 1))

.esp_empty <- function(nw, base, p, tabs) {
  if (!(base == -1L && p == 1L)) return(NULL)
  el <- .active_el(nw)
  if (nrow(el) == 0) return(NULL)
  A <- .active_adj(nw, el)
  sp <- (A %*% A)[el]
  colSums(tabs[sp + 1, , drop = FALSE])
}

.esp_signed_term <- function(nw, base, path, tabs, coef.names) {
  base <- .parse_base(base)
  p <- .parse_path(path)
  list(name = "esp_signed", coef.names = coef.names,
       inputs = as.double(tabs), iinputs = c(base, p),
       auxiliaries = ~.signaux, dependence = TRUE,
       emptynwstats = .esp_empty(nw, base, p, tabs))
}

.gwesp_signed_term <- function(nw, arglist, base, path, label, ...) {
  a <- check.ErgmTerm(nw, arglist, directed = FALSE, bipartite = FALSE,
                      varnames = c("decay"), vartypes = c("numeric"),
                      defaultvalues = list(NULL), required = c(TRUE))
  decay <- a$decay[1]
  n <- network.size(nw)
  tabs <- matrix(.gw_table(decay, max(n - 1, 1)), ncol = 1)
  .esp_signed_term(nw, base, path, tabs, paste0(label, ".fixed.", decay))
}

#' @templateVar name .signaux
#' @title Signed-network state (auxiliary)
#' @description Internal auxiliary that maintains the sign of every dyad and
#'   caches of shared partners by sign for the signed-network terms of
#'   \pkg{Bstergm}. It is requested automatically and never needs to be
#'   specified by the user.
#' @usage
#' # binary: .signaux
#' @name .signaux-ergmAuxiliary
#' @aliases InitErgmTerm..signaux
#' @keywords internal
#' @concept undirected
#' @concept signed
NULL
InitErgmTerm..signaux <- function(nw, arglist, ...) {
  check.ErgmTerm(nw, arglist, directed = FALSE, bipartite = FALSE)
  el <- .active_el(nw)
  list(name = "_signaux", coef.names = character(0),
       iinputs = c(nrow(el), el[, 1], el[, 2]), dependence = TRUE)
}

#' @name gwesp_signed-ergmTerm
#' @title Geometrically weighted edgewise shared partners by sign
#' @description Geometrically weighted edgewise shared partner statistic of a
#'   signed network, in which the base dyads have sign `base` and only shared
#'   partners whose two legs have the signs given by `path` are counted.
#'
#'   For a base dyad \eqn{\{i,j\}} let \eqn{sp_{ij}} be the number of nodes
#'   \eqn{k} such that the signs of \eqn{\{i,k\}} and \eqn{\{k,j\}} match
#'   `path`. The statistic is
#'   \deqn{e^{\alpha} \sum_{\{i,j\}:\, y_{ij} = base} \{1 - (1 - e^{-\alpha})^{sp_{ij}}\}}{exp(alpha) sum_{ij: y_ij = base} {1 - (1 - exp(-alpha))^sp_ij}}
#'
#'   It is evaluated on a signed network built by [as_signed_network()]: the
#'   ties of the network are the positive ties and the network attribute
#'   `.active` holds the interacting dyads, so that negative ties are active
#'   dyads that are not ties. The terms are implemented in C and are much
#'   faster than the equivalent multilayer
#'   [gwespL][ergm.multi::gwespL-ergmTerm] terms of \pkg{ergm.multi}.
#'
#' @usage
#' # binary: gwesp_signed(decay, base = "+", path = "++")
#' @param decay the decay parameter \eqn{\alpha} (required; the curved version
#'   is not available).
#' @param base sign of the base dyads, `"+"` or `"-"`.
#' @param path signs of the two legs of a shared partner path: `"++"`
#'   (shared friends), `"--"` (shared enemies) or `"+-"` (one positive and one
#'   negative leg).
#'
#' @details The special cases used by Caimo and Gollini are available as
#'   [gwesf_pos][gwesf_pos-ergmTerm] (\eqn{\mathrm{gwesf}^+}),
#'   [gwesf_neg][gwesf_neg-ergmTerm] (\eqn{\mathrm{gwesf}^-}),
#'   [gwese_pos][gwese_pos-ergmTerm] (\eqn{\mathrm{gwese}^+}) and
#'   [gwese_neg][gwese_neg-ergmTerm] (\eqn{\mathrm{gwese}^-}).
#'   The statistic \code{gwesp_signed(decay, "+", "++")} is identical to the
#'   ordinary \code{gwesp(decay, fixed = TRUE)} on the positive ties.
#' @seealso [as_signed_network()], [esp_signed][esp_signed-ergmTerm],
#'   [triangle_signed][triangle_signed-ergmTerm]
#' @examples
#' data(senators70_am)
#' nw <- as_signed_network(senators70_am$S99)
#' # The four balance-related statistics of Caimo and Gollini, and a mixed one
#' summary(nw ~ gwesf_pos(0.6) + gwese_pos(0.6) + gwesf_neg(0.6) + gwese_neg(0.6) +
#'           gwesp_signed(0.6, base = "+", path = "+-"))
#' # Simulate signs on the observed interactions
#' sim <- simulate(nw ~ edges + gwesf_pos(0.6) + gwese_neg(0.6), coef = c(-1, 0.3, -0.1),
#'                 constraints = ~active_signs, seed = 1)
#' table(signed_adjacency(sim)[upper.tri(diag(70))])
#' @aliases InitErgmTerm.gwesp_signed
#' @concept undirected
#' @concept signed
#' @concept triad-related
NULL
InitErgmTerm.gwesp_signed <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed = FALSE, bipartite = FALSE,
                      varnames = c("decay", "base", "path"),
                      vartypes = c("numeric", "character", "character"),
                      defaultvalues = list(NULL, "+", "++"),
                      required = c(TRUE, FALSE, FALSE))
  decay <- a$decay[1]
  n <- network.size(nw)
  tabs <- matrix(.gw_table(decay, max(n - 1, 1)), ncol = 1)
  base <- .parse_base(a$base)
  p <- .parse_path(a$path)
  lab <- paste0("gwesp_signed(", c("-", "", "+")[base + 2], ",", c("++", "--", "+-")[p + 1], ")")
  .esp_signed_term(nw, base, a$path, tabs, paste0(lab, ".fixed.", decay))
}

#' @name gwesf_pos-ergmTerm
#' @title Geometrically weighted positive edgewise shared friends
#' @description \eqn{\mathrm{gwesf}^+}: positive ties weighted by their number of
#'   shared friends (both legs positive). This is the balanced
#'   \eqn{(+,+,+)} triad viewed from a positive tie, and equals
#'   `gwesp(decay, fixed = TRUE)` on the positive ties. Equivalent to
#'   [gwesp_signed][gwesp_signed-ergmTerm]`(decay, base = "+", path = "++")`.
#' @usage
#' # binary: gwesf_pos(decay)
#' @param decay the decay parameter (required).
#' @seealso [gwesp_signed][gwesp_signed-ergmTerm], [as_signed_network()]
#' @aliases InitErgmTerm.gwesf_pos
#' @concept undirected
#' @concept signed
#' @concept triad-related
NULL
InitErgmTerm.gwesf_pos <- function(nw, arglist, ...)
  .gwesp_signed_term(nw, arglist, "+", "++", "gwesf_pos")

#' @name gwesf_neg-ergmTerm
#' @title Geometrically weighted negative edgewise shared friends
#' @description \eqn{\mathrm{gwesf}^-}: negative ties weighted by their number of
#'   shared friends (both legs positive), i.e. the unbalanced
#'   \eqn{(+,+,-)} triad viewed from its negative tie. Equivalent to
#'   [gwesp_signed][gwesp_signed-ergmTerm]`(decay, base = "-", path = "++")`.
#' @usage
#' # binary: gwesf_neg(decay)
#' @param decay the decay parameter (required).
#' @seealso [gwesp_signed][gwesp_signed-ergmTerm], [as_signed_network()]
#' @aliases InitErgmTerm.gwesf_neg
#' @concept undirected
#' @concept signed
#' @concept triad-related
NULL
InitErgmTerm.gwesf_neg <- function(nw, arglist, ...)
  .gwesp_signed_term(nw, arglist, "-", "++", "gwesf_neg")

#' @name gwese_pos-ergmTerm
#' @title Geometrically weighted positive edgewise shared enemies
#' @description \eqn{\mathrm{gwese}^+}: positive ties weighted by their number of
#'   shared enemies (both legs negative), i.e. the balanced \eqn{(+,-,-)}
#'   triad viewed from its positive tie. Equivalent to
#'   [gwesp_signed][gwesp_signed-ergmTerm]`(decay, base = "+", path = "--")` and
#'   to `gwespL(decay, fixed = TRUE, L.base = ~positive,
#'   Ls.path = list(~active & !positive, ~active & !positive))` in the
#'   multilayer representation.
#' @usage
#' # binary: gwese_pos(decay)
#' @param decay the decay parameter (required).
#' @seealso [gwesp_signed][gwesp_signed-ergmTerm], [as_signed_network()]
#' @aliases InitErgmTerm.gwese_pos
#' @concept undirected
#' @concept signed
#' @concept triad-related
NULL
InitErgmTerm.gwese_pos <- function(nw, arglist, ...)
  .gwesp_signed_term(nw, arglist, "+", "--", "gwese_pos")

#' @name gwese_neg-ergmTerm
#' @title Geometrically weighted negative edgewise shared enemies
#' @description \eqn{\mathrm{gwese}^-}: negative ties weighted by their number of
#'   shared enemies (both legs negative), i.e. the \eqn{(-,-,-)} triad, which
#'   is balanced under weak balance theory. Equivalent to
#'   [gwesp_signed][gwesp_signed-ergmTerm]`(decay, base = "-", path = "--")`.
#' @usage
#' # binary: gwese_neg(decay)
#' @param decay the decay parameter (required).
#' @seealso [gwesp_signed][gwesp_signed-ergmTerm], [as_signed_network()]
#' @aliases InitErgmTerm.gwese_neg
#' @concept undirected
#' @concept signed
#' @concept triad-related
NULL
InitErgmTerm.gwese_neg <- function(nw, arglist, ...)
  .gwesp_signed_term(nw, arglist, "-", "--", "gwese_neg")

#' @name esp_signed-ergmTerm
#' @title Edgewise shared partner distribution by sign
#' @description For each value in `d`, the number of dyads of sign `base` that
#'   have exactly that many shared partners whose legs match `path`. Useful
#'   mainly for goodness-of-fit checks.
#' @usage
#' # binary: esp_signed(d, base = "+", path = "++")
#' @param d a vector of non-negative integers.
#' @param base,path as in [gwesp_signed][gwesp_signed-ergmTerm].
#' @seealso [gwesp_signed][gwesp_signed-ergmTerm]
#' @aliases InitErgmTerm.esp_signed
#' @concept undirected
#' @concept signed
#' @concept triad-related
NULL
InitErgmTerm.esp_signed <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed = FALSE, bipartite = FALSE,
                      varnames = c("d", "base", "path"),
                      vartypes = c("numeric", "character", "character"),
                      defaultvalues = list(NULL, "+", "++"),
                      required = c(TRUE, FALSE, FALSE))
  d <- as.integer(a$d)
  if (any(d < 0)) stop("Argument 'd' must be non-negative.", call. = FALSE)
  n <- network.size(nw)
  L <- max(n - 1, 1)
  tabs <- vapply(d, function(dd) as.double((seq_len(L) - 1) == dd), numeric(L))
  tabs <- matrix(tabs, nrow = L)
  base <- .parse_base(a$base)
  p <- .parse_path(a$path)
  lab <- paste0("esp_signed(", c("-", "", "+")[base + 2], ",", c("++", "--", "+-")[p + 1], ")#", d)
  .esp_signed_term(nw, base, a$path, tabs, lab)
}

.degree_neg_term <- function(nw, tabs, coef.names) {
  el <- .active_el(nw)
  n <- network.size(nw)
  deg <- tabulate(c(el[, 1], el[, 2]), nbins = n)
  list(name = "degree_neg_signed", coef.names = coef.names,
       inputs = as.double(tabs), auxiliaries = ~.signaux, dependence = TRUE,
       emptynwstats = colSums(tabs[deg + 1, , drop = FALSE]))
}

#' @name gwdegree_neg-ergmTerm
#' @title Geometrically weighted degree of negative ties
#' @description The geometrically weighted degree statistic computed on the
#'   negative ties (active dyads that are not positive),
#'   \eqn{e^{\alpha}\sum_i \{1 - (1 - e^{-\alpha})^{d^-_i}\}}, where
#'   \eqn{d^-_i} is the number of negative ties of node \eqn{i}. For the
#'   positive ties use the ordinary `gwdegree(decay, fixed = TRUE)`.
#' @usage
#' # binary: gwdegree_neg(decay)
#' @param decay the decay parameter (required).
#' @seealso [degree_neg][degree_neg-ergmTerm], [as_signed_network()]
#' @aliases InitErgmTerm.gwdegree_neg
#' @concept undirected
#' @concept signed
NULL
InitErgmTerm.gwdegree_neg <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed = FALSE, bipartite = FALSE,
                      varnames = c("decay"), vartypes = c("numeric"),
                      defaultvalues = list(NULL), required = c(TRUE))
  decay <- a$decay[1]
  n <- network.size(nw)
  tabs <- matrix(.gw_table(decay, n), ncol = 1)
  .degree_neg_term(nw, tabs, paste0("gwdegree_neg.fixed.", decay))
}

#' @name degree_neg-ergmTerm
#' @title Degree distribution of negative ties
#' @description For each value in `d`, the number of nodes with exactly that
#'   many negative ties.
#' @usage
#' # binary: degree_neg(d)
#' @param d a vector of non-negative integers.
#' @seealso [gwdegree_neg][gwdegree_neg-ergmTerm]
#' @aliases InitErgmTerm.degree_neg
#' @concept undirected
#' @concept signed
NULL
InitErgmTerm.degree_neg <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed = FALSE, bipartite = FALSE,
                      varnames = c("d"), vartypes = c("numeric"),
                      defaultvalues = list(NULL), required = c(TRUE))
  d <- as.integer(a$d)
  if (any(d < 0)) stop("Argument 'd' must be non-negative.", call. = FALSE)
  n <- network.size(nw)
  tabs <- matrix(vapply(d, function(dd) as.double((seq_len(n) - 1) == dd), numeric(n)), nrow = n)
  .degree_neg_term(nw, tabs, paste0("degree_neg", d))
}

#' @name triangle_signed-ergmTerm
#' @title Signed triangle census
#' @description Number of triangles (three mutually active dyads) by sign
#'   composition. Under structural balance theory `"+++"` and `"+--"` are
#'   balanced, `"++-"` is unbalanced and `"---"` is balanced only under weak
#'   balance. Because the active dyads are fixed, the four counts sum to a
#'   constant, so at most three of them can enter a model together.
#' @usage
#' # binary: triangle_signed(type = c("+++", "++-", "+--", "---"))
#' @param type a character vector of sign patterns; the order of the signs
#'   within a pattern is irrelevant.
#' @seealso [gwesp_signed][gwesp_signed-ergmTerm], [signed_triad_census()]
#' @aliases InitErgmTerm.triangle_signed
#' @concept undirected
#' @concept signed
#' @concept triad-related
NULL
InitErgmTerm.triangle_signed <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed = FALSE, bipartite = FALSE,
                      varnames = c("type"), vartypes = c("character"),
                      defaultvalues = list(c("+++", "++-", "+--", "---")),
                      required = c(FALSE))
  type <- gsub("\\s", "", a$type)
  if (!all(grepl("^[+-]{3}$", type)))
    stop("Each 'type' must consist of three '+' or '-' signs, e.g. \"++-\".", call. = FALSE)
  npos <- nchar(gsub("-", "", type))
  canon <- vapply(npos, function(k) paste0(strrep("+", k), strrep("-", 3 - k)), "")
  empty <- NULL
  if (any(npos == 0)) {
    el <- .active_el(nw)
    ntri <- if (nrow(el)) {
      A <- .active_adj(nw, el)
      sum((A %*% A)[el]) / 3
    } else 0
    empty <- ifelse(npos == 0, ntri, 0)
  }
  list(name = "triangle_signed", coef.names = paste0("triangle_signed", canon),
       iinputs = as.integer(npos), auxiliaries = ~.signaux, dependence = TRUE,
       emptynwstats = empty)
}

#' @name active_signs-ergmConstraint
#' @title Only the signs of active dyads may change
#' @description Sample space constraint for signed networks built with
#'   [as_signed_network()]: only the active dyads (network attribute
#'   `.active`) may change, and only between positive and negative. This is
#'   the sample space of the sign process of a separable signed ERGM, where
#'   the interaction layer is held fixed.
#' @usage
#' # active_signs
#' @seealso [as_signed_network()]
#' @aliases InitErgmConstraint.active_signs
#' @concept signed
NULL
InitErgmConstraint.active_signs <- function(nw, arglist, ...) {
  check.ErgmTerm(nw, arglist)
  el <- .active_el(nw)
  list(free_dyads = function()
         as.rlebdm(as.edgelist(el, n = network.size(nw), directed = FALSE)),
       dependence = FALSE)
}
