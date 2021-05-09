#' calc_brier
#'
#' Simple brier score calculation
#'
#' @param predict \eqn{P(Y = 1)} for each value
#' @param outcome result (in \eqn{\{0, 1\}})
#' @param avg Whether to use the mean.
#' @return Length one double: the total brier sum, or the avg (default)
#' @export
calc_brier <- function(predict,
                       outcome,
                       avg = TRUE) {
  if (!avg) {
    return(sum((outcome - predict)^2))
  } else {
    return(mean((outcome - predict)^2))
  }
}
