#' Reshapes a list of simulations to a nice dataframe, by p-cut
#'
#' @inheritParams reshape_list_of_sims
#' @param p_cut Vector of one-sided permution p-values to use.
#' @return Data frame spreading results...
#'
#' @export
reshape_p_cut_list <- function(list_of_sims,
                               treat_model_name,
                               mu_model_name,
                               p_cut,
                               n_rows,
                               n_cols,
                               num_weight_vectors) {
  list_of_sims <- lapply(1L:length(list_of_sims), function(j) {
    if (is.null(list_of_sims[[j]][["id"]])) {
      list_of_sims[[j]][["id"]] <- j
    }
    return(list_of_sims[[j]])
  })

  do.call(rbind, lapply(p_cut, function(p_cut_val) {
    do.call(rbind, lapply(list_of_sims, function(par_res) {
      p_briers <- unlist(lapply(
        par_res[["weighted_results"]], function(x) {
          x[["permutation_brier"]]
        }
      ))
      if (all(p_briers > p_cut_val)) {
        given_cut_ind <- which.min(p_briers)
      } else {
        given_cut_ind <- which(p_briers <= p_cut_val)[1]
      }

      data.frame(
        id = par_res[["id"]],
        treat_model = treat_model_name,
        mu_model = mu_model_name,
        p_cut = p_cut_val,
        p_brier =
          par_res[["weighted_results"]][[given_cut_ind]][["permutation_brier"]],
        raw_brier =
          par_res[["weighted_results"]][[given_cut_ind]][["raw_brier"]],
        n_rows = n_rows,
        n_cols = n_cols,
        num_weight_vectors = num_weight_vectors,
        n_sinks = par_res[["weighted_results"]][[given_cut_ind]][["n_sinks"]],
        naive_est = par_res[["naive_est"]],
        propensity_est =
          par_res[["propensity_results"]][[given_cut_ind]][["est"]],
        mahal_est = par_res[["mahal_results"]][[given_cut_ind]][["est"]],
        weighted_est = par_res[["weighted_results"]][[given_cut_ind]][["est"]]
      )
    }))
  }))
}
