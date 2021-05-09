#' propensity_bipartite_matches
#'
#' Propensity match for bipartite matching.
#'
#' @inheritParams all_bipartite_matches
#' @param propensity_list See \code{gen_propensity_list}
#'
#' @export
propensity_bipartite_matches <- function(x_mat,
                                         treat_vec,
                                         match_method = c(
                                           "with_replacement",
                                           "optimal",
                                           "greedy"
                                         ),
                                         propensity_list =
                                           gen_propensity_list(),
                                         n_sinks = 0,
                                         caliper_list = gen_caliper_list(),
                                         sqrt_mahal = TRUE,
                                         tol_val = NULL) {
  ## in case of logical
  treat_vec <- treat_vec * 1L

  ## generate propensity score
  prop_score <- propensity_score(
    x_mat = x_mat,
    treat_vec = treat_vec,
    propensity_list = propensity_list
  )
  prop_dist_mat <- abs(outer(
    prop_score[treat_vec == 1],
    prop_score[treat_vec == 0],
    "-"
  ))

  if (!is.null(caliper_list)) {
    prop_dist_mat <- prop_dist_mat + create_caliper(caliper_list,
                                                    treat_vec = treat_vec
    )
  }

  bipartite_matches(
    dist_mat = prop_dist_mat,
    treat_vec = treat_vec,
    match_method = match_method,
    n_sinks = n_sinks,
    tol_val = tol_val
  )
}
