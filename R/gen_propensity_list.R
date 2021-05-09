#' gen_propensity_list
#'
#' Takes in elements needed for propensity function, checks input and builds a named list.
#'
#' @param propensity_function A function that accepts a list with four elements: \code{x_train}, \code{x_test}, \code{y_train}, \code{y_test}, and forms predictions using \code{x_test} (I guess \code{y_test} isn't used yet)
#' @param oos_propensity Logical, whether to predict out of sample for the propensity score.
#' @param n_folds Default NULL; how many folds if using out of sample propensity.
#' @return Named list, same names as input params.
#'
#' @export
gen_propensity_list <- function(propensity_function = propensity_score_xgb(),
                                oos_propensity = FALSE,
                                n_folds = NULL) {
  stopifnot(is_tf(oos_propensity))

  if (missing(oos_propensity) && !is.null(n_folds)) {
    oos_propensity <- TRUE
  }

  if (!oos_propensity) {
    if (!is.null(n_folds)) {
      stop("n_folds shouldn't be set if not using out of sample")
    }
  } else {
    n_folds <- ifelse(is.null(n_folds), 5L, n_folds)
    stopifnot(is.numeric(n_folds) && length(n_folds) == 1L &&
                !is.na(n_folds))
  }

  list(
    propensity_function = propensity_function,
    oos_propensity = oos_propensity,
    n_folds = n_folds
  )
}
