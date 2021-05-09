#' compute_sim_result
#'
#' Takes in functions to generate simulation data, and computes simulation results for our method vs. Mahalanobis and propensity score matching.
#'
#' @inheritParams generate_simulation_input
#' @inheritParams all_bipartite_matches
#' @param n_sink_gen Default \code{n_sink_generator}. This argument should be a function that accepts a \code{treat_vec} and produces a vector of sink numbers.
#' @param num_weight_vectors How many weight vectors to generate.
#' @param silent Whether to suppress messages as it's running.. Default \code{!interactive()}.
#' @return Returns a named list:
#' \describe{
#'     \item{\code{naive_est}}{Mean difference between all treated units and all control}
#'     \item{\code{propensity_results}}{List of lists: each with \code{n_sinks} and the \code{est}}
#'     \item{\code{mahal_results}}{Same as above}
#'     \item{\code{weighted_results}}{List of lists: each with \code{n_sinks}, the raw brier score, the permutation brier score, and the \code{est}}
#' }
#'
#' @export
compute_sim_result <- function(x_generator = default_x_generator,
                               treat_prob_generator,
                               mean_generator,
                               error_generator = default_error_generator,
                               n_sink_gen = n_sink_generator(),
                               match_method = "with_replacement",
                               n_rows = 500L,
                               n_cols = 5L,
                               num_weight_vectors = 100L,
                               silent = !interactive()) {
  sim_data <- generate_simulation_input(
    n_rows = n_rows,
    n_cols = n_cols,
    x_generator = x_generator,
    treat_prob_generator = treat_prob_generator,
    mean_generator = mean_generator,
    error_generator = error_generator
  )

  x_mat <- sim_data[["x_mat"]]
  y_vector <- sim_data[["y_vec"]]
  treat_vec <- sim_data[["treat_vec"]]

  rm(sim_data)

  n_sinks <- n_sink_gen(treat_vec)
  match_list_est_func <- (function(y_vector, treat_vec) {
    function(match_list) {
      match_estimate(
        match_list = match_list,
        y_vector = y_vector,
        treat_vec = treat_vec
      )
    }
  })(y_vector, treat_vec)

  list_est_func <- (function(n_sinks) {
    function(match_lists) {
      Map(
        function(n_sink, match_list) {
          list(
            n_sinks = n_sink,
            est = match_list_est_func(match_list)
          )
        },
        n_sinks,
        match_lists
      )
    }
  })(n_sinks)

  naive_est <- mean(y_vector[treat_vec == 1]) - mean(y_vector[treat_vec == 0])

  ## ------------------------------------

  if (!silent) {
    message("propensity matches")
  }

  propensity_matches <- propensity_bipartite_matches(
    x_mat = x_mat,
    treat_vec = treat_vec,
    match_method = match_method,
    propensity_list = gen_propensity_list(
      propensity_function = propensity_score_linear(),
      oos_propensity = FALSE
    ),
    n_sinks = n_sinks
  )

  propensity_ests <- list_est_func(propensity_matches)

  ## ------------------------------------

  if (!silent) {
    message("mahal matches")
  }

  mahal_matches <- all_bipartite_matches(
    x_mat = x_mat,
    cov_x = covariance_with_ranks(x_mat),
    weight_list = list(rep(1 / n_cols, times = n_cols)),
    treat_vec = treat_vec,
    match_method = match_method,
    n_sinks = n_sinks
  )[[1]]

  mahal_ests <- list_est_func(mahal_matches)

  ## ------------------------------------

  if (!silent) {
    message("all weighted matches")
  }

  weight_list <- generate_random_weights(
    prior_weights = rep(1 / n_cols, times = n_cols),
    number_vectors = num_weight_vectors,
    minimum_weights = rep(1 / (3 * n_cols), times = n_cols)
  )

  brier_matches <- brier_bipartite_matches(
    x_mat = x_mat,
    cov_x = covariance_with_ranks(x_mat),
    weight_list = weight_list,
    treat_vec = treat_vec,
    match_method = match_method,
    n_sinks = n_sinks,
    silent = silent
  )

  permutation_results <- permutation_matches(
    matches_by_sinks = brier_matches[["matches_by_sinks"]],
    briers_by_sinks = brier_matches[["briers_by_sinks"]],
    x_mat = x_mat,
    n_sinks = n_sinks,
    silent = silent
  )

  weighted_results <- lapply(
    permutation_results[["best_matches"]],
    function(match_results) {
      list(
        n_sinks = match_results[["n_sinks"]],
        raw_brier = match_results[["raw_brier"]],
        permutation_brier = match_results[["permutation_brier"]],
        est = match_list_est_func(match_results[["match_list"]])
      )
    }
  )

  ## ------------------------------------

  list(
    naive_est = naive_est,
    propensity_results = propensity_ests,
    mahal_results = mahal_ests,
    weighted_results = weighted_results
  )
}
