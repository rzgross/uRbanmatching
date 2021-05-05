#' exp / 1 + exp
#' @param x numeric vector
#' @return numeric vector, now in (0, 1)
expit <- function(x) {
  exp(x) / (1 + exp(x))
}
