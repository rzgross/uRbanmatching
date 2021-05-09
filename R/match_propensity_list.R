#' match_propensity_list
#'
#' Generates the propensity parameters used for using propensity-based calipers.
#'
#' @inheritParams gen_propensity_list
#' @param caliper_sd_mult Maximum gap between units.In \code{sd(propensity_score) * k}, where this parameter is the value k. Default 0.6.
#' @param continuous_mult Penality for going above a certain caliper.
#' @return list with names equal to all input params
#' @export
match_propensity_list <- function(propensity_function = propensity_score_xgb(),
                                  oos_propensity = FALSE,
                                  n_folds = NULL,
                                  caliper_sd_mult = 0.6,
                                  continuous_mult = 100) {
  if (is.null(propensity_function)) {
    return(NULL)
  }

  plain_prop_list <- gen_propensity_list(
    propensity_function = propensity_function,
    oos_propensity = oos_propensity,
    n_folds = n_folds
  )

  list(
    propensity_function = plain_prop_list[["propensity_function"]],
    oos_propensity = plain_prop_list[["oos_propensity"]],
    n_folds = plain_prop_list[["n_folds"]],
    caliper_sd_mult = caliper_sd_mult,
    continuous_mult = continuous_mult
  )
}
