#' paper_treatment_functions
#'
#' Generates the four different treatment functions used for simulations in the paper.
#'
#' @param target_mean Desired mean for each of the treatment probabilities, should be in (0, 1).
#' @return List of four functions that take in a matrix and output a vector of probabilities:
#' \describe{
#'   \item{\code{constant_treat_prob}}{Just repeats \code{target_mean}}
#'   \item{\code{logistic_treat_prob}}{Computes a logistic style relationship
#'   between the data matrix and the output vector}
#'   \item{\code{sparse_treat_prob}}{Based on the sign of the first data matrix
#'   column; gives a number a bit above the mean and a bit below.}
#'   \item{\code{sparse_nonlin_treat_prob}}{Cubic logit function of the
#'   first vector}
#' }
#'
#' @export
paper_treatment_functions <- function(target_mean = 0.425) {
  stopifnot(length(target_mean) == 1 &&
              0 < target_mean && target_mean < 1)

  constant_treat_prob <- function(x_mat) {
    rep(target_mean, nrow(x_mat))
  }

  logistic_treat_prob <- function(x_mat) {
    coef_vec <- rnorm(ncol(x_mat), 0, 0.2)
    lin_vec <- c(x_mat %*% coef_vec)

    target_mean_expit(
      target_mean = target_mean,
      linear_vector = lin_vec
    )
  }

  sparse_treat_prob <- function(x_mat) {
    median_adjusted <- x_mat[, 1] - median(x_mat[, 1])

    min_gap <- min(target_mean, 1 - target_mean) / 2

    ifelse(sign(median_adjusted) > 0,
           target_mean + min_gap,
           target_mean - min_gap
    )
  }

  sparse_nonlin_treat_prob <- function(x_mat) {
    mean_adj <- (x_mat[, 1] - mean(x_mat[, 1])) / (sd(x_mat[, 1]) * 2) + 1
    numer <- mean_adj^3 - mean_adj^2

    target_mean_expit(
      target_mean = target_mean,
      linear_vector = numer
    )
  }

  ## ------------------------------------

  list(
    constant_treat_prob = constant_treat_prob,
    logistic_treat_prob = logistic_treat_prob,
    sparse_treat_prob = sparse_treat_prob,
    sparse_nonlin_treat_prob = sparse_nonlin_treat_prob
  )
}
