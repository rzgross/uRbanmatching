#' weighted_mahal
#'
#' Computes weighted Mahalanobis distance, using Choleski decomp.
#'
#' @param x_mat Numeric matrix (adjust non-numeric columns prior), already rank-adjusted if desired.
#' @param cov_x Covariance of x, calculated potentially with ranks.
#' @param weight_vec Vector of weights corresponding to columns of x_mat, giving weights relative to normal Mahalanobis.
#' @param treat_vec optionally specify which units are treated.
#' @param sqrt_mahal Logical, default TRUE; whether to use regular Mahalanobis.
#' @param partial_index How to partition the full matrix.
#' @return returns a matrix of pairwise distances; the relevant indexing depends on \code{treat_vec} and \code{partial_index}
#' @export
weighted_mahal <- function(x_mat,
                           cov_x,
                           weight_vec = NULL,
                           treat_vec = NULL,
                           sqrt_mahal = TRUE,
                           partial_index = NULL) {
  if (!is.null(partial_index) && !is.null(treat_vec)) {
    stop("Supply at most one of treat_vec and partial_index", call. = FALSE)
  }
  if (!is.matrix(x_mat)) {
    stop("x_mat must be a matrix", call. = FALSE)
  }

  chol_cov <- chol(cov_x)

  if (is.null(weight_vec)) {
    weight_vec <- rep(1 / ncol(x_mat), times = ncol(x_mat))
  } else {
    if (length(weight_vec) != ncol(x_mat)) {
      stop("`weight_vec` should have length equal to `ncol(x_mat)`")
    }
  }

  weighted_ymat <- forwardsolve(diag(1 / weight_vec) %*%
                                  t(chol_cov), t(x_mat))

  ymat_sumsq <- colSums(weighted_ymat^2)

  if (!is.null(treat_vec) || !is.null(partial_index)) {
    if (!is.null(treat_vec)) {
      ## in case it's logical...:
      treat_vec <- treat_vec * 1

      row_index <- treat_vec == 1
      col_index <- treat_vec == 0
    } else {
      row_index <- partial_index[[1]]
      col_index <- partial_index[[2]]
    }

    crossprod_weighted <- t(weighted_ymat[, row_index]) %*%
      weighted_ymat[, col_index]

    mahal_mat <- outer(ymat_sumsq[row_index],
                       ymat_sumsq[col_index],
                       FUN = "+"
    ) -
      2 * crossprod_weighted
  } else {
    crossprod_weighted <- t(weighted_ymat) %*% weighted_ymat

    mahal_mat <- outer(ymat_sumsq, ymat_sumsq, FUN = "+") -
      2 * crossprod_weighted
  }

  ## correcting for e.g. system tol (well, working with doubles)
  mahal_mat <- pmax(mahal_mat, 0)

  if (sqrt_mahal) {
    mahal_mat <- sqrt(mahal_mat)
  }

  mahal_mat
}
