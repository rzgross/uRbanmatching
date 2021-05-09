#' optimal_match
#'
#' Computes optimal matches. Alternative to "greedy" and "with replacement".
#'
#' @inheritParams bipartite_matches
#' @param n_sinks single value
#' @export
optimal_match <- function(dist_mat,
                          treat_vec,
                          n_sinks = 0,
                          tol_val = 1e-4) {
  ## ------------------------------------

  dist_mat_sink <- cbind(
    dist_mat,
    ## sinks are zeros
    matrix(0, nrow = nrow(dist_mat), ncol = n_sinks)
  )
  main_treat_index <- which(treat_vec == 1L)
  rownames(dist_mat_sink) <- main_treat_index

  all_treated <- treat_vec
  if (n_sinks > 0) {
    ## the sinks are not treated
    all_treated <- c(treat_vec, rep(0L, n_sinks))
  }
  all_control_index <- which(all_treated == 0L)
  colnames(dist_mat_sink) <- all_control_index

  treat_data <- data.frame(treated = all_treated)
  rownames(treat_data) <- seq_len(nrow(treat_data))

  ## Matching using optmatch
  ## if need be:
  ## options("optmatch_max_problem_size" = Inf))
  match_vec <- optmatch::pairmatch(dist_mat_sink,
                                   tol = tol_val,
                                   data = treat_data
  )
  list_results <- aggregate(names(match_vec),
                            by = list(match_vec),
                            FUN = function(x) {
                              as.numeric(x)
                            },
                            simplify = FALSE
  )[["x"]]
  first_vec <- unlist(lapply(list_results, function(x) x[1]))
  second_vec <- unlist(lapply(list_results, function(x) x[2]))
  control_ind <- all_treated[first_vec] == 0
  keep_ind <- pmax(first_vec, second_vec) <= length(treat_vec)

  treat_subj_ind <- ifelse(control_ind, second_vec, first_vec)[keep_ind]
  control_subj_ind <- ifelse(control_ind, first_vec, second_vec)[keep_ind]

  ord_vec <- order(treat_subj_ind)

  treat_subj_ind <- treat_subj_ind[ord_vec]
  control_subj_ind <- control_subj_ind[ord_vec]

  match_list <- list(
    treat_index = treat_subj_ind,
    treat_index_within = match(treat_subj_ind, main_treat_index),
    control_index = control_subj_ind,
    control_index_within = match(control_subj_ind, all_control_index)
  )
  match_list[["distance"]] <- dist_mat[cbind(
    match_list[["treat_index_within"]],
    match_list[["control_index_within"]]
  )]

  match_list
}
