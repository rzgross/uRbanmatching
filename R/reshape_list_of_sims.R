#' Reshapes a list of simulations to a nice dataframe
#'
#' In the future this should be cleaner
#'
#' @param list_of_sims List of results from \code{compute_sim_results}.
#' @param treat_model_name Name of the treatment model.
#' @param mu_model_name Name of the mean generation model.
#' @param n_rows How many rows were used.
#' @param n_cols How many columns were used.
#' @param num_weight_vectors How many weight vectors were used
#' @return Data frame spreading results...
#' @author Colman Humphrey
#'
#' @export
reshape_list_of_sims <- function(list_of_sims,
                                 treat_model_name,
                                 mu_model_name,
                                 n_rows,
                                 n_cols,
                                 num_weight_vectors) {
  list_of_sims <- lapply(1L:length(list_of_sims), function(j) {
    if (is.null(list_of_sims[[j]][["id"]])) {
      list_of_sims[[j]][["id"]] <- j
    }
    return(list_of_sims[[j]])
  })

  do.call(rbind, lapply(list_of_sims, function(sim_res) {
    do.call(rbind, lapply(1L:length(sim_res[["weighted_results"]]), function(j) {
      data.frame(
        id = sim_res[["id"]],
        treat_model = treat_model_name,
        mu_model = mu_model_name,
        p_brier =
          sim_res[["weighted_results"]][[j]][["permutation_brier"]],
        raw_brier =
          sim_res[["weighted_results"]][[j]][["raw_brier"]],
        n_rows = n_rows,
        n_cols = n_cols,
        num_weight_vectors = num_weight_vectors,
        n_sinks = sim_res[["weighted_results"]][[j]][["n_sinks"]],
        naive_est = sim_res[["naive_est"]],
        propensity_est =
          sim_res[["propensity_results"]][[j]][["est"]],
        mahal_est = sim_res[["mahal_results"]][[j]][["est"]],
        weighted_est = sim_res[["weighted_results"]][[j]][["est"]]
      )
    }))
  }))
}
