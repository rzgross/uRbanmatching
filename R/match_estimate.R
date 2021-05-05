#' Computes a simple mean difference in an outcome vector
#' between treatment and control in a paired match.
#' Computes the average difference between the treated units and the
#' control units for a match, given as a match list.  Optionally
#' confirms that all treated units are indeed treated.
#' @param match_list Typical \code{match_list} object from
#'     \code{bipartite_matches}.
#' @param y_vector The outcome vector.
#' @param treat_vec Default NULL, provide if you want it checked.
#' @return Returns a single number, the mean difference.
#'
#' @export
match_estimate <- function(match_list,
                           y_vector,
                           treat_vec = NULL) {
  if (!is.null(treat_vec)) {
    stopifnot(length(treat_vec) == length(y_vector))
    stopifnot(all(treat_vec[match_list[["treat_index"]]] == 1L))
    stopifnot(all(treat_vec[match_list[["control_index"]]] == 0L))
  }

  mean(y_vector[match_list[["treat_index"]]] -
         y_vector[match_list[["control_index"]]])
}
