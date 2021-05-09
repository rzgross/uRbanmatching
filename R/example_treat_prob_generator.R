#' example_treat_prob_generator
#'
#' Function that generates treatment probabilities
#' given an input matrix.
#'
#' @param x_mat Numeric matrix.
#' @return A vector of probabilities, with length equal to \code{nrow(x_mat)}.
#'
#' @export
example_treat_prob_generator <- function(x_mat) {
  flat_1 <- (x_mat[, 1] - mean(x_mat[, 1])) / sd(x_mat[, 1])
  numer <- (flat_1^3 - 2 * flat_1^2) / 10

  if (ncol(x_mat) > 1) {
    numer <- numer + sign(x_mat[, 2]) / 2 +
      x_mat %*% rnorm(
        n = ncol(x_mat),
        mean = 0, sd = 0.1
      )
  }

  expit(numer)
}
