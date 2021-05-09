#' paper_mean_functions
#'
#' Generates the four different mean functions used for simulations in the paper.
#'
#' @return List of four functions that take in a matrix and
#'   output a vector of means:
#' \describe{
#'   \item{\code{constant_mu}}{Vector of zeroes}
#'   \item{\code{linear_mu}}{The matrix multiplied by a random vector.}
#'   \item{\code{sign_mu}}{Literally \code{sign(x_mat[, 1])}}
#'   \item{\code{non_linear_mu}}{Messy non-linear function of the data matrix}
#' }
#'
#' @export
paper_mean_functions <- function() {
  constant_mu <- function(x_mat) {
    rep(0, nrow(x_mat))
  }

  linear_mu <- function(x_mat) {
    coef_vec <- rnorm(ncol(x_mat), 0, 1)
    lin_vec <- c(x_mat %*% coef_vec)

    lin_vec - mean(lin_vec)
  }

  sign_mu <- function(x_mat) {
    median_adjusted <- x_mat[, 1] - median(x_mat[, 1])

    sign(median_adjusted) - mean(sign(median_adjusted))
  }

  non_linear_mu <- function(x_mat) {
    coef_vec <- rnorm(ncol(x_mat), 0, 0.5)
    min_abs <- apply(abs(x_mat), 1, min)

    end_col <- ncol(x_mat)

    lin_vec <- cos(x_mat[, 1]) +
      sin(c(x_mat %*% coef_vec)) * sign(x_mat[, end_col]) -
      pmin(tan(pi / 3 - min_abs * pi), 10)

    lin_vec - mean(lin_vec)
  }

  ## ------------------------------------

  list(
    constant_mu = constant_mu,
    linear_mu = linear_mu,
    sign_mu = sign_mu,
    non_linear_mu = non_linear_mu
  )
}
