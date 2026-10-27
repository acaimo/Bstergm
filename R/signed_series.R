#' Sequence of signed networks
#'
#' Creates the data object used by [bstergm_sign()] (and [bstergm()]): a sequence of undirected
#' signed networks \eqn{y_0, \dots, y_T} observed on the same set of nodes.
#'
#' @param nets a list of at least two signed networks, each given as a signed
#'   adjacency matrix (entries -1, 0, 1), a `network` with a `sign` edge
#'   attribute, a network built by [as_signed_network()] or a two-layer
#'   [ergm.multi::Layer()] network with layers `positive` and `active`.
#' @param vertex.attr optional data frame (or named list) of vertex attributes,
#'   one row per node, made available to the model terms (e.g.
#'   `nodematch("Party")`).
#' @param sign.attr name of the edge attribute holding the signs when the
#'   networks are `network` objects.
#'
#' @details Each transition \eqn{y_{t-1} \to y_t} is decomposed, for both the
#'   interaction layer \eqn{x} (\eqn{x_{ij} = |y_{ij}|}) and the sign layer
#'   \eqn{z} (\eqn{z_{ij} = 1} if \eqn{y_{ij} = +1}), into
#'   \describe{
#'     \item{formation}{\eqn{x^F_t = x_{t-1} \cup x_t}, \eqn{z^F_t = z_{t-1} \cup z_t};}
#'     \item{persistence}{\eqn{x^P_t = x_{t-1} \cap x_t}, \eqn{z^P_t = z_{t-1} \cap z_t}.}
#'   }
#'   These networks are returned by [transition_networks()].
#'
#' @return An object of class `signed_series`: a list with the signed
#'   adjacency matrices `y`, the number of nodes `n`, the time labels `times`
#'   and the vertex attributes `vertex.attr`.
#' @seealso [bstergm_sign()], [bstergm()], [transition_networks()]
#' @examples
#' data(senators70_am)
#' senate <- signed_series(senators70_am[1:3],
#'                         vertex.attr = attr(senators70_am$S99, "senators"))
#' senate
#' summary(senate)
#' @export
signed_series <- function(nets, vertex.attr = NULL, sign.attr = "sign") {
  if (inherits(nets, "signed_series")) {
    if (!is.null(vertex.attr)) nets$vertex.attr <- .check_vertex_attr(vertex.attr, nets$n)
    return(nets)
  }
  if (is.matrix(nets) || is.network(nets))
    stop("'nets' must be a list of at least two signed networks.", call. = FALSE)
  nets <- as.list(nets)
  if (length(nets) < 2) stop("At least two time points are needed.", call. = FALSE)

  if (is.null(vertex.attr)) vertex.attr <- .vertex_attr_from(nets[[1]])
  Y <- lapply(nets, signed_adjacency, sign.attr = sign.attr)
  n <- unique(vapply(Y, nrow, 1L))
  if (length(n) != 1) stop("All networks must have the same number of nodes.", call. = FALSE)

  node.names <- rownames(Y[[1]])
  if (is.null(node.names)) node.names <- as.character(seq_len(n))
  Y <- lapply(Y, function(m) {
    dimnames(m) <- list(node.names, node.names)
    m
  })
  times <- names(nets)
  if (is.null(times) || any(!nzchar(times))) times <- paste0("t", seq_along(nets) - 1)
  names(Y) <- times

  structure(list(y = Y, n = n, times = times, node.names = node.names,
                 vertex.attr = .check_vertex_attr(vertex.attr, n)),
            class = "signed_series")
}

.vertex_attr_from <- function(x) {
  if (!is.network(x)) return(NULL)
  if (.is_layer_network(x)) x <- uncombine_network(x)[[1]]
  va <- setdiff(list.vertex.attributes(x), c("na", "vertex.names", ".LayerID", ".LayerName"))
  va <- va[!startsWith(va, ".")]
  if (!length(va)) return(NULL)
  as.data.frame(lapply(stats::setNames(va, va), function(a) get.vertex.attribute(x, a)),
                stringsAsFactors = FALSE)
}

.check_vertex_attr <- function(vertex.attr, n) {
  if (is.null(vertex.attr)) return(NULL)
  va <- as.data.frame(vertex.attr, stringsAsFactors = FALSE)
  if (nrow(va) != n) stop("'vertex.attr' must have one row per node (", n, ").", call. = FALSE)
  rownames(va) <- NULL
  va
}

#' @rdname signed_series
#' @param x a `signed_series` object.
#' @param ... unused.
#' @export
print.signed_series <- function(x, ...) {
  if (isTRUE(x$binary)) {
    cat("Sequence of", length(x$y), "undirected binary networks on", x$n, "nodes\n")
    tab <- .series_table(x)[, c("time", "active", "density")]
    names(tab)[2] <- "ties"
    print(tab, row.names = FALSE)
  } else {
    cat("Sequence of", length(x$y), "undirected signed networks on", x$n, "nodes\n")
    print(.series_table(x), row.names = FALSE)
  }
  if (!is.null(x$vertex.attr))
    cat("Vertex attributes:", paste(names(x$vertex.attr), collapse = ", "), "\n")
  invisible(x)
}

.series_table <- function(x) {
  nd <- x$n * (x$n - 1) / 2
  data.frame(time = x$times,
             positive = vapply(x$y, function(m) sum(m > 0) / 2, 0),
             negative = vapply(x$y, function(m) sum(m < 0) / 2, 0),
             active = vapply(x$y, function(m) sum(m != 0) / 2, 0),
             density = round(vapply(x$y, function(m) sum(m != 0) / 2, 0) / nd, 3),
             `%positive` = round(100 * vapply(x$y, function(m) sum(m > 0) / max(sum(m != 0), 1), 0), 1),
             check.names = FALSE)
}

#' @rdname signed_series
#' @param object a `signed_series` object.
#' @return `summary()` returns a list with tables of tie counts by time point,
#'   tie counts of the formation and persistence networks of each transition,
#'   and the signed triad census by time point.
#' @export
summary.signed_series <- function(object, ...) {
  tn <- transition_networks(object)
  trans <- data.frame(
    transition = paste(head(object$times, -1), tail(object$times, -1), sep = " -> "),
    F.active = vapply(tn$F, function(m) sum(m != 0) / 2, 0),
    F.positive = vapply(tn$F, function(m) sum(m > 0) / 2, 0),
    P.active = vapply(tn$P, function(m) sum(m != 0) / 2, 0),
    P.positive = vapply(tn$P, function(m) sum(m > 0) / 2, 0),
    check.names = FALSE)
  out <- list(times = .series_table(object), transitions = trans,
              triads = t(vapply(object$y, signed_triad_census, integer(4))),
              binary = isTRUE(object$binary))
  if (out$binary) {
    out$times <- out$times[, c("time", "active", "density")]
    names(out$times)[2] <- "ties"
    out$transitions <- out$transitions[, c("transition", "F.active", "P.active")]
    names(out$transitions)[2:3] <- c("F.ties", "P.ties")
    out$triads <- out$triads[, "+++", drop = FALSE]
    colnames(out$triads) <- "triangles"
  }
  class(out) <- "summary.signed_series"
  out
}

#' @export
print.summary.signed_series <- function(x, ...) {
  cat("Ties by time point:\n")
  print(x$times, row.names = FALSE)
  cat("\nFormation (F) and persistence (P) networks by transition:\n")
  print(x$transitions, row.names = FALSE)
  cat(if (isTRUE(x$binary)) "\nTriangles by time point:\n" else "\nSigned triangles by time point:\n")
  print(x$triads)
  invisible(x)
}

#' Formation and persistence networks of each transition
#'
#' @param data a [signed_series()] object.
#' @return A list with components `F` and `P`, each a list of signed adjacency
#'   matrices (one per transition), and `prev`, the signed adjacency matrices
#'   of the previous time points.
#' @seealso [signed_series()]
#' @export
transition_networks <- function(data) {
  data <- signed_series(data)
  Y <- data$y
  Tn <- length(Y) - 1
  form <- pers <- vector("list", Tn)
  for (t in seq_len(Tn)) {
    a <- Y[[t]]; b <- Y[[t + 1]]
    xF <- (a != 0) | (b != 0); zF <- (a > 0) | (b > 0)
    xP <- (a != 0) & (b != 0); zP <- (a > 0) & (b > 0)
    form[[t]] <- .signed_from_layers(zF & xF, xF, dimnames(a))
    pers[[t]] <- .signed_from_layers(zP & xP, xP, dimnames(a))
  }
  nm <- paste(head(data$times, -1), tail(data$times, -1), sep = "->")
  names(form) <- names(pers) <- nm
  list(F = form, P = pers, prev = Y[seq_len(Tn)], cur = Y[seq_len(Tn) + 1])
}

.signed_from_layers <- function(z, x, dn = NULL) {
  m <- matrix(0L, nrow(x), ncol(x), dimnames = dn)
  m[x] <- -1L
  m[z] <- 1L
  m
}

#' @rdname signed_series
#' @param layout optional matrix of node coordinates; by default a common
#'   layout of the union of all interactions is used.
#' @param vertex.col vertex colours, or the name of a vertex attribute to colour
#'   by (levels in alphabetical order get grey shades from dark to light).
#' @param edge.col grey levels of positive (solid) and negative (dotted) ties.
#' @export
plot.signed_series <- function(x, layout = NULL, vertex.col = NULL,
                               edge.col = c(positive = "grey15", negative = "grey60"), ...) {
  K <- length(x$y)
  any_act <- Reduce(`|`, lapply(x$y, function(m) m != 0)) * 1
  if (is.null(layout)) {
    g <- network::as.network(any_act, directed = FALSE)
    layout <- network::network.layout.fruchtermanreingold(g, NULL)
  }
  if (is.character(vertex.col) && length(vertex.col) == 1 && !is.null(x$vertex.attr) &&
      vertex.col %in% names(x$vertex.attr)) {
    f <- factor(x$vertex.attr[[vertex.col]])
    vertex.col <- grDevices::gray.colors(nlevels(f), start = 0.1, end = 0.9)[as.integer(f)]
  }
  if (is.null(vertex.col)) vertex.col <- "grey30"
  edge.col <- rep_len(edge.col, 2)
  op <- par(mfrow = c(1, K), mar = c(0.5, 0.5, 2, 0.5))
  on.exit(par(op))
  for (k in seq_len(K)) {
    m <- x$y[[k]]
    plot(layout, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, main = x$times[k], ...)
    e <- which(upper.tri(m) & m != 0, arr.ind = TRUE)
    if (nrow(e)) {
      s <- m[e]
      graphics::segments(layout[e[, 1], 1], layout[e[, 1], 2], layout[e[, 2], 1], layout[e[, 2], 2],
                         col = ifelse(s > 0, adjustcolor(edge.col[1], 0.45), adjustcolor(edge.col[2], 0.25)),
                         lty = ifelse(s > 0, 1, 3))
    }
    points(layout, pch = 21, bg = vertex.col, col = "black", cex = 1.3)
  }
  invisible(layout)
}

#' Sequence of binary networks
#'
#' Creates the data object used by [bstergm()]: a sequence of undirected
#' binary networks \eqn{y_0, \dots, y_T} observed on the same set of nodes.
#'
#' @param nets a list of at least two undirected binary networks, given as
#'   adjacency matrices or `network` objects. A [signed_series()] object is
#'   returned unchanged: its interactions (ties of either sign) are then the
#'   binary networks.
#' @param vertex.attr optional data frame (or named list) of vertex attributes,
#'   one row per node; by default taken from the first network if it is a
#'   `network` object.
#'
#' @return An object of class `c("network_series", "signed_series")` with
#'   all ties coded as 1.
#' @seealso [bstergm()], [signed_series()]
#' @examples
#' data(senators70_am)
#' # Interactions (ties of either sign) of the Senate data as binary networks
#' nets <- lapply(senators70_am[1:3], function(m) abs(m))
#' network_series(nets, vertex.attr = attr(senators70_am$S99, "senators"))
#' @export
network_series <- function(nets, vertex.attr = NULL) {
  if (inherits(nets, "signed_series")) {
    if (!is.null(vertex.attr)) nets$vertex.attr <- .check_vertex_attr(vertex.attr, nets$n)
    return(nets)
  }
  if (is.matrix(nets) || is.network(nets))
    stop("'nets' must be a list of at least two networks.", call. = FALSE)
  nets <- as.list(nets)
  if (is.null(vertex.attr)) vertex.attr <- .vertex_attr_from(nets[[1]])
  Y <- lapply(nets, function(x) {
    if (is.network(x)) {
      if (is.directed(x)) stop("Only undirected networks are supported.", call. = FALSE)
      m <- as.matrix(x, matrix.type = "adjacency")
      dimnames(m) <- list(network.vertex.names(x), network.vertex.names(x))
    } else m <- as.matrix(x)
    if (anyNA(m)) stop("Missing values are not supported.", call. = FALSE)
    (m != 0) * 1L
  })
  out <- signed_series(Y, vertex.attr = vertex.attr)
  out$binary <- TRUE
  class(out) <- c("network_series", class(out))
  out
}
