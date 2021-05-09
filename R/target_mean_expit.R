#' target_mean_expit
#'
#' Shifts the linear component until the mean of the expit is the target mean.
#'
#' @param target_mean Desired target mean.
#' @param linear_vector Basically \eqn{X \beta} (no shift).
#' @return \code{expit(linear_vector + a)} for some \code{a} such that \code{mean(expit(linear_vector + a))} is close to \code{target_mean}.
#'
#' @export
target_mean_expit <- function(target_mean,
                              linear_vector) {
  mono_func <- function(alpha_value) {
    mean(expit(linear_vector + alpha_value))
  }

  alpha_value <- binary_search(
    target_value = target_mean,
    monotone_function = mono_func
  )

  expit(linear_vector + alpha_value)
}
