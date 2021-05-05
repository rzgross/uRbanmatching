#' Converts \code{rank_cols} in all allowed forms to an integer index
#'
#' Takes in ranked cols, either given by integer (numeric is fine) index;
#' logical index; or named. Converts to integer index
#' @param rank_cols integer/number, or logical, or names within
#'   \code{colnames(x_mat)}
#' @param x_mat the x matrix of interest
#' @keywords internal
rank_integer_index <- function(rank_cols, x_mat) {
  if (is.null(rank_cols)) {
    return(vector("integer", 0L))
  }

  if (all(rank_cols %in% 1L:ncol(x_mat))) {
    return(rank_cols)
  }

  if (is.logical(rank_cols)) {
    if (length(rank_cols) == ncol(x_mat)) {
      return(which(rank_cols))
    } else {
      stop("logical `rank_cols` must be same length as `ncol(x_mat)`")
    }
  } else {
    if (is.null(colnames(x_mat))) {
      stop("x_mat must have colnames to use named rank_cols")
    }

    if (any(!(rank_cols %in% colnames(x_mat)))) {
      stop("not all rank_cols are present in x_mat colnames")
    }

    return(which(colnames(x_mat) %in% rank_cols))
  }
}
