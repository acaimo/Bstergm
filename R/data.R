#' Signed networks of 70 U.S. Senators, 99th to 102nd Congress
#'
#' Signed networks among the 70 senators who served continuously in the 99th
#' (1985-1987), 100th (1987-1989), 101st (1989-1991) and 102nd (1991-1993)
#' U.S. Congresses. Positive and negative ties were extracted from bill
#' co-sponsorship with the stochastic degree sequence model backbone of Neal
#' (2014, 2020): a positive tie means that two senators co-sponsored
#' significantly more bills than expected, a negative tie significantly fewer.
#'
#' @format A list of 4 signed adjacency matrices (`S99`, `S100`, `S101`,
#'   `S102`), each 70 x 70 and symmetric with entries -1, 0 and 1. Each matrix
#'   carries the attribute `senators`, a data frame with columns `ID`,
#'   `Senator` (name), `State` and `Party` (`"D"` or `"R"`), which can be
#'   passed as the `vertex.attr` argument of [signed_series()].
#'
#' @source Neal, Z. P. (2020). A sign of the times? Weak and strong
#'   polarization in the U.S. Congress, 1973-2016. \emph{Social Networks}, 60,
#'   103-112. Data:
#'   <https://figshare.com/articles/dataset/A_Sign_of_the_Times/8096429?file=20758701>
#' @references Neal, Z. P. (2014). The backbone of bipartite projections:
#'   Inferring relationships from co-authorship, co-sponsorship, co-attendance
#'   and other co-behaviors. \emph{Social Networks}, 39, 84-97.
#' @examples
#' data(senators70_am)
#' sapply(senators70_am, function(m) c(positive = sum(m > 0) / 2, negative = sum(m < 0) / 2))
#' table(attr(senators70_am$S99, "senators")$Party)
"senators70_am"
