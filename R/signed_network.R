#' Signed network representation used by the native signed terms
#'
#' Converts a signed network to the representation used by the signed-network
#' ergm terms of this package: a [network][network::network] whose ties are
#' the *positive* ties, with the active (interacting) dyads stored as a
#' two-column edge list in the network attribute `.active`. Negative ties are
#' the active dyads that are not ties.
#'
#' @param y a signed network: a square symmetric matrix with entries in
#'   \{-1, 0, 1\}; an undirected [network][network::network] whose edges carry
#'   the edge attribute `sign` (or a different attribute given by `sign.attr`);
#'   or a two-layer [ergm.multi::Layer()] network with layers `positive` and
#'   `active` as used by earlier versions of this package.
#' @param vertex.attr optional data frame (or named list) of vertex attributes
#'   with one row per node.
#' @param sign.attr name of the edge attribute holding the signs when `y` is a
#'   network.
#' @param ... unused.
#'
#' @return An undirected `network` of positive ties with network attribute
#'   `.active`.
#' @seealso [signed_adjacency()], [signed_triad_census()],
#'   [gwese_pos][gwese_pos-ergmTerm]
#' @examples
#' data(senators70_am)
#' y <- senators70_am$S99
#' nw <- as_signed_network(y, vertex.attr = attr(y, "senators"))
#' nw
#' summary(nw ~ edges + gwesf_pos(0.6) + gwese_pos(0.6) + gwese_neg(0.6) +
#'           triangle_signed())
#' @export
as_signed_network <- function(y, vertex.attr = NULL, sign.attr = "sign", ...) {
  m <- signed_adjacency(y, sign.attr = sign.attr)
  if (is.null(vertex.attr) && is.network(y) && !.is_layer_network(y)) {
    va <- list.vertex.attributes(y)
    va <- setdiff(va, c("na", "vertex.names"))
    if (length(va)) vertex.attr <- lapply(stats::setNames(va, va), function(a) get.vertex.attribute(y, a))
  }
  .signed_network_from_matrix(m, vertex.attr)
}

.signed_network_from_matrix <- function(m, vertex.attr = NULL, pos = NULL) {
  n <- nrow(m)
  ut <- upper.tri(m)
  act <- which(ut & m != 0, arr.ind = TRUE)
  if (is.null(pos)) pos <- which(ut & m > 0, arr.ind = TRUE)
  nw <- .network_from_el(pos, n)
  storage.mode(act) <- "integer"
  dimnames(act) <- NULL
  nw %n% ".active" <- act
  if (!is.null(rownames(m))) network.vertex.names(nw) <- rownames(m)
  .set_vertex_attr(nw, vertex.attr)
}

.network_from_el <- function(el, n) {
  nw <- network.initialize(n, directed = FALSE)
  el <- as.matrix(el)
  if (length(el) && nrow(el)) {
    nw <- network.edgelist(cbind(as.integer(el[, 1]), as.integer(el[, 2])), nw)
  }
  nw
}

.set_vertex_attr <- function(nw, vertex.attr) {
  if (is.null(vertex.attr)) return(nw)
  vertex.attr <- as.list(vertex.attr)
  n <- network.size(nw)
  for (a in names(vertex.attr)) {
    v <- vertex.attr[[a]]
    if (is.factor(v)) v <- as.character(v)
    if (length(v) != n) stop("Vertex attribute ", sQuote(a), " must have one value per node.", call. = FALSE)
    nw <- set.vertex.attribute(nw, a, v)
  }
  nw
}

.is_layer_network <- function(y) {
  is.network(y) && inherits(y, "combined_networks") && ".LayerName" %in% list.vertex.attributes(y)
}

#' Signed adjacency matrix
#'
#' Returns the signed adjacency matrix (entries -1, 0, 1) of a signed network
#' in any of the forms accepted by [as_signed_network()].
#'
#' @inheritParams as_signed_network
#' @return A square symmetric integer matrix.
#' @export
signed_adjacency <- function(y, sign.attr = "sign") {
  if (is.network(y)) {
    if (is.directed(y)) stop("Only undirected signed networks are supported.", call. = FALSE)
    n <- network.size(y)
    act <- y %n% ".active"
    if (!is.null(act)) {
      m <- matrix(0L, n, n)
      act <- as.matrix(act)
      if (length(act)) m[act] <- m[act[, 2:1, drop = FALSE]] <- -1L
      pos <- as.matrix(y, matrix.type = "edgelist")
      if (length(pos)) m[pos[, 1:2, drop = FALSE]] <- m[pos[, 2:1, drop = FALSE]] <- 1L
      dimnames(m) <- list(network.vertex.names(y), network.vertex.names(y))
      return(m)
    }
    if (.is_layer_network(y)) {
      layers <- uncombine_network(y)
      nm <- names(layers)
      if (!all(c("positive", "active") %in% nm))
        stop("A multilayer network must have layers named 'positive' and 'active'.", call. = FALSE)
      P <- as.matrix(layers[[match("positive", nm)]])
      A <- as.matrix(layers[[match("active", nm)]])
      m <- A * (2L * P - 1L)
      m[P == 1 & A == 0] <- 1L
      storage.mode(m) <- "integer"
      return(m)
    }
    if (!(sign.attr %in% list.edge.attributes(y)))
      stop("The network has no edge attribute ", sQuote(sign.attr), " with the tie signs.", call. = FALSE)
    m <- as.matrix(y, matrix.type = "adjacency", attrname = sign.attr)
    return(.check_signed_matrix(m))
  }
  .check_signed_matrix(y)
}

.check_signed_matrix <- function(m) {
  if (inherits(m, "Matrix")) m <- as.matrix(m)
  m <- as.matrix(m)
  if (nrow(m) != ncol(m)) stop("A signed adjacency matrix must be square.", call. = FALSE)
  if (anyNA(m)) stop("Missing values in signed networks are not supported.", call. = FALSE)
  if (!all(m %in% c(-1, 0, 1))) stop("Signed adjacency matrices must only contain -1, 0 and 1.", call. = FALSE)
  if (!isSymmetric(unname(m) * 1)) stop("Only undirected (symmetric) signed networks are supported.", call. = FALSE)
  diag(m) <- 0
  out <- matrix(as.integer(m), nrow(m), ncol(m), dimnames = dimnames(m))
  out
}

#' Signed triad census
#'
#' Counts the triangles of a signed network by their number of positive ties.
#' Under structural balance theory, triangles with three or one positive ties
#' are balanced; under weak balance, all-negative triangles are balanced too.
#'
#' @inheritParams as_signed_network
#' @return A named integer vector with the counts of `"+++"`, `"++-"`,
#'   `"+--"` and `"---"` triangles.
#' @examples
#' data(senators70_am)
#' sapply(senators70_am, signed_triad_census)
#' @export
signed_triad_census <- function(y, sign.attr = "sign") {
  m <- signed_adjacency(y, sign.attr = sign.attr)
  P <- (m > 0) * 1
  N <- (m < 0) * 1
  tr <- function(A, B, C) sum(diag(A %*% B %*% C))
  out <- c("+++" = tr(P, P, P) / 6, "++-" = tr(P, P, N) / 2,
           "+--" = tr(P, N, N) / 2, "---" = tr(N, N, N) / 6)
  storage.mode(out) <- "integer"
  out
}
