#' binary_search
#'
#' Generic binary search algo: finding the input.
#'
#' @param target_value The value the function should achieve.
#' @param monotone_function The function we want to find the  relevant input for.
#' @param init_bounds Default NULL; Bounds on the function being searched.
#' @param error_gap Binary search until the gap between the function on the input and the target value is smaller than this number.
#' @param max_iters Maximum umber of iterations to try.
#' @return The input that gives \code{target_value} as output.
#' @export

binary_search <- function(target_value,
                          monotone_function,
                          init_bounds = NULL,
                          error_gap = 1e-6,
                          max_iters = 100L) {
  stopifnot(error_gap > 0)

  test_bounds <- init_bounds
  if (is.null(init_bounds)) {
    test_bounds <- c(0, 1)
  }
  if (monotone_function(test_bounds[1L]) >
      monotone_function(test_bounds[2L])) {
    return(binary_search(
      target_value = -target_value,
      monotone_function = function(x) {
        -monotone_function(x)
      },
      init_bounds = init_bounds,
      max_iters = max_iters
    ))
  }

  if (is.null(init_bounds)) {
    if (monotone_function(0) > target_value) {
      upper_bound <- 0
      lower_bound <- -1
      while (monotone_function(lower_bound) > target_value) {
        lower_bound <- lower_bound * 1.1 - 1
      }
    } else {
      lower_bound <- 0
      upper_bound <- 1
      while (monotone_function(upper_bound) < target_value) {
        upper_bound <- upper_bound * 1.1 + 1
      }
    }
  } else {
    stopifnot(length(init_bounds) == 2L)
    stopifnot(init_bounds[1L] < init_bounds[2L])

    lower_bound <- init_bounds[1]
    upper_bound <- init_bounds[2]
  }

  stopifnot(sign(monotone_function(lower_bound) - target_value) *
              sign(monotone_function(upper_bound) - target_value) == -1L)

  input_val <- (lower_bound + upper_bound) / 2
  mono_val <- monotone_function(input_val)

  iters <- 0L

  while (abs(mono_val - target_value) > error_gap && iters < max_iters) {
    if (mono_val > target_value) {
      upper_bound <- input_val
    } else {
      lower_bound <- input_val
    }

    input_val <- (lower_bound + upper_bound) / 2
    mono_val <- monotone_function(input_val)

    iters <- iters + 1L
  }

  if (iters == max_iters) {
    stop("reached max iterations")
  }

  input_val
}
