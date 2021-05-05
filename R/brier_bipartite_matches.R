#' Computes all matches, then gets the brier scores for each. Reorder by
#' number of sinks.
#'
#' @inheritParams all_bipartite_matches
#' @inheritParams brier_score_cv
#' @return List of matches within sink values,
#'  and brier scores for each.
#' @param match_method ***
#' @param propensity_list ***
#' @param sqrt_mahal ***
#' @param silent ***
#'
#' @export
brier_bipartite_matches <- function(x_mat,
                                    cov_x,
                                    weight_list,
                                    treat_vec,
                                    match_method = c(
                                      "with_replacement",
                                      "optimal",
                                      "greedy"
                                    ),
                                    n_sinks = 0L,
                                    caliper_list = gen_caliper_list(),
                                    propensity_list =
                                      match_propensity_list(NULL),
                                    sqrt_mahal = TRUE,
                                    tol_val = NULL,
                                    design = "cross_all",
                                    num_folds = 5,
                                    match_predict_function =
                                      match_predict_xgb(),
                                    silent = !interactive()) {
  if (is.null(n_sinks)) {
    n_sinks <- 0L
  }

  ## generate all matches: one per weight vector per n_sink value
  all_matches <- all_bipartite_matches(
    x_mat = x_mat,
    cov_x = cov_x,
    weight_list = weight_list,
    treat_vec = treat_vec,
    match_method = match_method,
    n_sinks = n_sinks,
    caliper_list = caliper_list,
    propensity_list = propensity_list,
    sqrt_mahal = sqrt_mahal,
    tol_val = tol_val
  )

  if (!silent) {
    message("getting briers")
  }

  ## get all brier scores for all results
  briers_by_sinks <- lapply(all_matches, function(all_by_sink) {
    if (!silent) {
      print(all_by_sink[[1]]["num_sinks"])
    }
    unlist(lapply(all_by_sink, function(indiv_match_list) {
      brier_score_cv(
        x_mat = x_mat,
        match_list = indiv_match_list,
        design = design,
        num_folds = num_folds,
        match_predict_function = match_predict_function
      )
    }))
  })

  list(
    matches_by_sinks = all_matches,
    briers_by_sinks = briers_by_sinks
  )
}
