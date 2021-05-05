#' Computes weighted Mahalanobis distance, using Choleski decomp.
#'
#' Note: def. of weighted:
#' \eqn{d(x_i, x_j) = (x_i - x_j)' W \Sigma^{-1} W (x_i - x_j)}
#' where \eqn{W} = \code{diag(weight_vec)}
#' R's cholesky gives \eqn{U} s.t. \eqn{U' U = S} (i.e. \eqn{U = chol(S)})
#' so in general we want:
#' \deqn{
#'   (x_i - x_j)' W (U' U)^{-1} W (x_i - x_j) \\
#' = (x_i - x_j)' W U^{-1} (U')^{-1} W (x_i - x_j) \\
#' = x_i' W U^{-1} (U')^{-1} W x_i +
#'   x_j' W U^{-1} (U')^{-1} W x_j \\
#'   - 2 x_i' W U^{-1} (U')^{-1} W x_j
#' }
#' Solving the above is manageable if we have \eqn{y_i = (U')^{-1} W x_i}
#' or moving terms around, \eqn{W^{-1} U' y_i = x_i},
#' which is simple by \code{forwardsolve}. Then we have:
#' \deqn{dist(x_i, x_j) = ||y_i||^2 + ||y_j||^2 - 2 y_i' y_j}
#' Letting \eqn{Y = (y_1', y_2', .... ,y_n')'}
#' (which we can get in one line,
#' \eqn{Y' = } \code{forwardsolve} \eqn{(W^{-1} U', x_mat')})
#' the first two parts are just \code{rowSums}\eqn{(Y^2)}
#' (note that the code uses \eqn{Y'}, thus \code{colSums})
#' (add \code{outer(.,.)} to finish)
#' and the last is \eqn{Y Y'}
#'
#' @param x_mat numeric matrix (adjust non-numeric columns prior),
#'   already rank-adjusted if desired
#' @param cov_x covariance of x, calculated potentially with ranks
#' @param weight_vec vector of weights corresponding to columns of x_mat,
#'   giving weights relative to "raw" (ranked) Mahalanobis.
#'   Note that the resulting matrix does depend on the scale of the
#'   weight vector, but the matching won't: scaling the Mahalanobis
#'   matrix has no effect on distance minimising pairs etc
#' @param treat_vec optionally specify which units are treated.
#'   If \code{NULL} (default), will just return
#'   \code{nrow(x_mat) x nrow(x_mat)} distance
#'   matrix of all pairs. Can be logicals or \eqn{{0, 1}}
#' @param sqrt_mahal logical, default TRUE; do you want regular Mahalanobis:
#'   \eqn{d(x_i, x_j) = (x_i - x_j)' \Sigma^{-1} (x_i - x_j)}
#'   or the square root? (in weighted, \eqn{\Sigma^{-1}} becomes
#'   \eqn{W \Sigma^{-1} W} in the above)
#' @param partial_index In some cases, you want a subset of the full
#'   N x N matrix, but you can't partition it like you want
#'   with treat_vec. e.g. you want
#'   \code{full_dist[c(1, 2, 3), 1:10]}
#'   then use
#'   \code{partial_index = list(c(1,2,3), 1:10)}
#' @return returns a matrix of pairwise distances; the relevant indexing
#'   depends on \code{treat_vec} and \code{partial_index}
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
