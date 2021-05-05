#' An example function that generates a mean vector
#' given an input matrix.
#'
#' @param x_mat Numeric matrix.
#' @return A vector of means, with length equal to \code{nrow(x_mat)}.
#'
#' @export
example_mean_generator <- function(x_mat) {
  b_vec <- rnorm(
    n = ncol(x_mat),
    mean = 0,
    sd = 0.3
  )

  linear_part <- x_mat %*% b_vec

  square_part <- (x_mat[, 1] - mean(x_mat[, 1]))^2
  cross_part <- 0
  if (ncol(x_mat) > 1) {
    square_part <- square_part - (x_mat[, 2] - mean(x_mat[, 2]))^2
    cross_part <- (x_mat[, 1] - mean(x_mat[, 1])) *
      (x_mat[, ncol(x_mat)] - mean(x_mat[, ncol(x_mat)]))
  }

  linear_part + square_part + cross_part
}
