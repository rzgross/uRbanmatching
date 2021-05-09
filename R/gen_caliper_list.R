#' gen_caliper_list
#'
#' Wrapper function to unify caliper input.
#'
#' @param caliper_vec Default NULL; numeric vector that stops matches if beyond \code{caliper_max}.
#' @param caliper_max The maximum allowed difference.
#' @param continuous_mult The value to multiply differences above caliper max. Set as \code{Inf} to have infinite penalties.
#' @return Either \code{NULL}, or a list with the same names as the input, after checking values.
#'
#' @export
gen_caliper_list <- function(caliper_vec = NULL,
                             caliper_max = NULL,
                             continuous_mult = 100) {
  if (is.null(caliper_vec)) {
    if (!is.null(caliper_max)) {
      stop("can't give `caliper_max` without `caliper_vec`")
    }
    return(NULL)
  }

  if (is.null(caliper_max)) {
    stop("supply `caliper_max` if using calipers")
  }

  if (length(caliper_max) > 1L) {
    stop("`caliper_max` should be length one")
  }

  if (length(continuous_mult) > 1L) {
    stop("`continuous_mult` should be length one")
  }

  return(list(
    caliper_vec = caliper_vec,
    caliper_max = caliper_max,
    continuous_mult = continuous_mult
  ))
}
