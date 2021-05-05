#' Converts indicated columns to ranked versions of themselves
#'
#' This function takes a numeric matrix and converts any columns
#' indicated by the \code{rank_cols} input to their ranks
#' (break ties however you want), scaled down by \code{nrow(x_mat)}.
#' @param x_mat numeric matrix (adjust non-numeric columns prior)
#' @param rank_cols names or index of columns to be converted to ranks
#'   before analysis
#' @param ties_method how to break ties in ranks, by default uses
#'   "average" (same as the rank function's default). See \code{?rank}
#'   for further options
#' @return x_mat again, potentially adjusted for ranks
#' @keywords internal
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
