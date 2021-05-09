#' optimal_sink_wrap
#'
#' Given a vector of sink values, generates an optimal match for each.
#'
#' @inheritParams bipartite_matches
#' @param n_sinks Default NULL, vector of sink values to use.
#' @param weight_vec Default \code{NULL}: optionally supply the weight vector used to generate \code{dist_mat} and it'll be returned in the \code{match_list} generated from this function.
#' @return List of lists.
#'
#' @export
optimal_sink_wrap <- function(dist_mat,
                              treat_vec,
                              n_sinks,
                              tol_val,
                              weight_vec = NULL) {
  if (is.null(n_sinks)) {
    n_sinks <- 0L
  }

  ## ------------------------------------

  setNames(lapply(n_sinks, function(sink_val) {
    match_list <- optimal_match(
      dist_mat,
      treat_vec,
      sink_val,
      tol_val
    )
    match_list[["num_sinks"]] <- sink_val
    if (!is.null(weight_vec)) {
      match_list[["weight_vec"]] <- weight_vec
    }
    match_list
  }), n_sinks)
}
