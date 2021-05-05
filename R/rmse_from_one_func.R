#' Calculates RMSE with a known target of one
#'
#' @param vec Vector of estimates around 1.
#' @return Returns a single number: the RMSE
#'
#' @export
rmse_from_one_func <- function(vec) {
  sqrt(mean((vec - 1)^2))
}
