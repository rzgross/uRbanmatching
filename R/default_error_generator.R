#' default_error_generator
#'
#' Function to generate mean-zero noise. Used as default in \code{generate_simulation_input}.
#'
#' @param n_rows Number of rows to produce.
#' @return Returns just mean 0 variance 1 normal noise.
#'
#' @export
default_error_generator <- function(n_rows) {
  rnorm(n = n_rows)
}
