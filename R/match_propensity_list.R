#' Generates the propensity parameters used for using propensity-based calipers
#'
#' We use this for input for \code{all_propensity_caliper_matches} (and
#' likely a nonbipartite version soon).
#' @inheritParams gen_propensity_list
#' @param caliper_sd_mult We'll set the maximum gap between units
#'   as \code{sd(propensity_score) * k}, where this parameter is the
#'   value k. Default 0.6.
#' @param continuous_mult See e.g. \code{gen_caliper_list}: instead of
#'   blocking matches that are "too far apart" on the caliper, we'll
#'   add a penalty for going above.
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
