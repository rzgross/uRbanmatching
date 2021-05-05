#' Wraps a match_list from a simple match result and slices off
#' various numbers of sinks
#'
#' Note that for \code{with_replacement} matching, this is the
#' right thing to do, but for greedy matching with a random ordering
#' it isn't obviously correct, but does give the nice property
#' of more sinks = lower distances.
#' @param simple_match_list match result from one of \code{with_replacement_match},
#'   \code{greedy_match}, \code{with_replacement_nbp_match} or
#'   \code{greedy_nbp_match}.
#' @param n_sinks default NULL, vector of sink values to use.
#' @param weight_vec Default \code{NULL}: optionally supply the weight vector
#'   used to generate \code{dist_mat} and it'll be returned in the
#'   \code{match_list} generated from this function
#' @return list of lists; see parent function
#'
#' @keywords internal
simple_sink_wrap <- function(simple_match_list,
                             n_sinks = NULL,
                             weight_vec = NULL) {
  if (is.null(n_sinks)) {
    n_sinks <- 0L
  }

  ## ------------------------------------

  dist_ranks <- rank(simple_match_list[["distance"]],
                     ties.method = "random"
  )
  setNames(lapply(n_sinks, function(sink_val) {
    keep_ind <- dist_ranks <= (length(dist_ranks) - sink_val)
    match_list <- lapply(simple_match_list, function(x) {
      x[keep_ind]
    })
    match_list[["num_sinks"]] <- sink_val
    if (!is.null(weight_vec)) {
      match_list[["weight_vec"]] <- weight_vec
    }
    match_list
  }), n_sinks)
}
