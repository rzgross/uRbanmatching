#' Default function to generate mean-zero noise. Very simple.
#'
#' Exported just so all functions have a default in
#' \code{generate_simulation_input}.
#' @param n_rows How many rows to produce.
#' @return Returns just mean 0 variance 1 normal noise.
#'
#' @export
default_error_generator <- function(n_rows) {
  rnorm(n = n_rows)
}
