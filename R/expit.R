#' expit
#'
#' Computes exp / 1 + exp, which is used in several functions.
#'
#' @param x numeric vector
#' @return numeric vector, now in (0, 1)
#' @export
expit <- function(x) {
  exp(x) / (1 + exp(x))
}
