#' all_bipartite_matches
#'
#' @param x_mat Input/design matrix
#' @param cov_x The covariance of \code{x_mat}. Can be rank_adjusted by using \code{covariance_with_ranks}.
#' @param weight_list List of weight vectors. See `generate_random_weights` to
#'   automatically generate a reasonable set of vectors.
#' @param treat_vec Logical (or 1/0) vector, indicating treatment (or control).
#' @param n_sinks Number of potential matches that don't need to be matched.
#' @param caliper_list Forces matches that are close on some metric.
#' @param tol_val For optimal matches, you can set a tolerance to be within optimality of, which can be zero for perfect optimality.
#' @param match_method "with_replacement", "optimal", or "greedy"
#' @param propensity_list Default NULL. List of propensity scores.
#' @param sqrt_mahal Whether to use square root of mahalanobis distances.
#' @import stats
#' @export

all_bipartite_matches <- function(x_mat,
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
                                  propensity_list = match_propensity_list(NULL),
                                  sqrt_mahal = TRUE,
                                  tol_val = NULL) {
  if (is.null(n_sinks)) {
    n_sinks <- 0L
  }

  if (!is.null(propensity_list)) {
    if (!is.null(caliper_list)) {
      stop(
        "don't use both `caliper_list` and `propensity_list`: ",
        " If you do want both, create the combined caliper separately"
      )
    }

    ## in case of logical
    treat_vec <- treat_vec * 1L

    ## generate propensity score
    prop_list_names <- c(
      "propensity_function",
      "oos_propensity",
      "n_folds"
    )
    prop_score <- propensity_score(
      x_mat = x_mat,
      treat_vec = treat_vec,
      propensity_list = propensity_list[prop_list_names]
    )
    caliper_list <- gen_caliper_list(
      caliper_vec = prop_score,
      caliper_max = sd(prop_score) * propensity_list[["caliper_sd_mult"]],
      continuous_mult = propensity_list[["continuous_mult"]]
    )
  }

  if (!is.null(caliper_list)) {
    caliper_dist_mat <- create_caliper(caliper_list,
                                       treat_vec = treat_vec
    )
  }

  by_weight_list <- lapply(weight_list, function(weight_vec) {
    w_dist_mat <- weighted_mahal(x_mat,
                                 cov_x = cov_x,
                                 weight_vec = weight_vec,
                                 treat_vec = treat_vec,
                                 sqrt_mahal = sqrt_mahal
    )

    if (!is.null(caliper_list)) {
      w_dist_mat <- w_dist_mat + caliper_dist_mat
    }

    bipartite_matches(
      dist_mat = w_dist_mat,
      treat_vec = treat_vec,
      match_method = match_method,
      n_sinks = n_sinks,
      tol_val = tol_val,
      weight_vec = weight_vec
    )
  })

  ## more natural to group by sink value
  setNames(lapply(n_sinks, function(x) {
    lapply(by_weight_list, function(y) {
      y[[as.character(x)]]
    })
  }), n_sinks)
}
