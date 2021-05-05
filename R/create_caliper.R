#' Creates caliper penalties
#'
#' Given a vector, \code{caliper_vec}, this function sets penalties
#' for any pairwise differences above \code{caliper_max}; either
#' an \code{Inf} penalty (the default), or a continuous penalty. Typical usage will
#' then add the resulting matrix from this function onto a distance
#' matrix, say from pairwise Mahalanobis.
#'
#' @param caliper_list Result of \code{gen_caliper_list}
#' @param treat_vec Optional; if you only want pairs between treat and control.
#' @return A matrix, either square (\code{|caliper_vec| x |caliper_vec|})
#'   or else \code{sum(treat_vec == 1) x sum(treat_vec == 0)}.
#'
#' @export
create_caliper <- function(caliper_list,
                           treat_vec = NULL) {
  if (is.null(treat_vec)) {
    caliper_treat <- caliper_list[["caliper_vec"]]
    caliper_control <- caliper_list[["caliper_vec"]]
  } else {
    if (!(length(treat_vec) == length(caliper_list[["caliper_vec"]]))) {
      stop("treat_vec not the same length as `caliper_vec`")
    }
    ## in case it's logical...:
    treat_vec <- treat_vec * 1L

    caliper_treat <- caliper_list[["caliper_vec"]][treat_vec == 1]
    caliper_control <- caliper_list[["caliper_vec"]][treat_vec == 0]
  }

  abs_caliper_diff <- abs(outer(caliper_treat,
                                caliper_control,
                                FUN = "-"
  )) - caliper_list[["caliper_max"]]
  abs_caliper_diff[abs_caliper_diff < 0] <- 0

  continuous_mult <- caliper_list[["continuous_mult"]]

  if (!is.null(continuous_mult) && continuous_mult < Inf) {
    penalty_mat <- continuous_mult * abs_caliper_diff
  } else {
    penalty_mat <- abs_caliper_diff
    penalty_mat[penalty_mat > 0] <- Inf
  }

  penalty_mat
}
