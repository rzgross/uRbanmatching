#' Simple brier score calc
#'
#' @param predict \eqn{P(Y = 1)} for each value
#' @param outcome result (in \eqn{\{0, 1\}})
#' @param avg logica, default TRUE, do you want the mean?
#' @return length one double: the total brier sum, or the avg (default)
#' @keywords internal
calc_brier <- function(predict,
                       outcome,
                       avg = TRUE) {
  if (!avg) {
    return(sum((outcome - predict)^2))
  } else {
    return(mean((outcome - predict)^2))
  }
}
