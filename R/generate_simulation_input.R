#' generate_simulation_input
#'
#' Generates input needed for a simulation run.
#'
#' @param n_rows How many rows to generate.
#' @param n_cols How many columns to use.
#' @param x_generator Function that takes number of rows and columns
#'   and produces a data matrix of that dimension.
#' @param treat_prob_generator Function that takes a matrix and produces
#'   treatment probabilities.
#' @param mean_generator Function that takes a matrix and produces
#'   an expected value for each row.
#' @param error_generator Function that accepts a number of rows and generates an error vector (e.g. normal noise).
#' @return List:
#' \describe{
#'   \item{\code{x_mat}}{Data Matrix}
#'   \item{\code{treat_vec}}{Treatment Vector}
#'   \item{\code{y_vec}}{Output Vector}
#' }
#'
#' @export
generate_simulation_input <- function(n_rows = 500L,
                                      n_cols = 5L,
                                      x_generator =
                                        default_x_generator,
                                      treat_prob_generator =
                                        example_treat_prob_generator,
                                      mean_generator =
                                        example_mean_generator,
                                      error_generator =
                                        default_error_generator) {
  x_mat <- x_generator(
    n_rows = n_rows,
    n_cols = n_cols
  )

  treat_prob <- treat_prob_generator(x_mat)
  treat_vec <- rbinom(length(treat_prob), 1, treat_prob)

  mu_pre_treat <- mean_generator(x_mat)

  list(
    x_mat = x_mat,
    treat_vec = treat_vec,
    y_vec = mu_pre_treat + treat_vec + error_generator(n_rows = n_rows)
  )
}
