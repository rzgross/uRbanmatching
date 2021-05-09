#' ranked_x
#'
#' Converts indicated columns to ranked versions of themselves. Takes a numeric matrix and converts any columns indicated by the \code{rank_cols} input to their ranks.
#'
#' @param x_mat Numeric matrix (adjust non-numeric columns prior)
#' @param rank_cols Names or index of columns to be converted to ranks
#' @param ties_method How to break ties in ranks, by default uses "average" (same as the rank function's default).
#'   for further options
#' @return x_mat, potentially adjusted for ranks
#' @export
ranked_x <- function(x_mat,
                     rank_cols = NULL,
                     ties_method = "average") {
  ## convert to ranks
  for (rank_change in rank_integer_index(rank_cols, x_mat)) {
    ## divides by nrow, better scaling
    x_mat[, rank_change] <- rank(x_mat[, rank_change],
                                 ties.method = ties_method
    ) / nrow(x_mat)
  }

  x_mat
}
