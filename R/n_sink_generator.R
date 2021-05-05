#' Generates a function that generates a vector of sink lengths
#'
#' This function returns a function that accepts a treatment vector
#' and generates a vector of numbers to be used as sink counts.
#' @param start_frac Smaller fraction of units to use as sink number, default 0.
#' @param end_frac Larger fraction of units to use as sink number, default 0.8.
#' @param length_out How many sink values we want, default 9.
#' @return Function that accepts a \code{treat_vec}
#'   and returns a vector of numbers.
#'
#' @export
n_sink_generator <- function(start_frac = 0,
                             end_frac = 0.8,
                             length_out = 9) {
  stopifnot(length(start_frac) == 1L &&
              length(end_frac) == 1L &&
              length(length_out) == 1L)
  stopifnot(0 <= start_frac && start_frac <= 1)
  stopifnot(0 <= end_frac && end_frac <= 1)
  stopifnot(length_out >= 1L)

  if (length_out == 1) {
    stopifnot(start_frac == end_frac)
  } else {
    stopifnot(start_frac < end_frac)
  }


  function(treat_vec) {
    num_treat <- sum(treat_vec)

    floor(seq(
      from = floor(start_frac * num_treat),
      to = floor(end_frac * num_treat),
      length.out = length_out
    ))
  }
}
