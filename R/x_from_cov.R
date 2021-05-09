#' x_from_cov
#'
#' Generates a matrix where each row is generated with covariance equal to \code{cov_mat}.
#'
#' @param cov_mat The covariance desired.
#' @param n_rows Number of rows to produce.
#' @return A matrix such that the covariance would tend to \code{cov_mat} as \code{n_rows} grows.
#'
#' @export
x_from_cov <- function(cov_mat,
                       n_rows) {
  ## recall: in R, chol gives U s.t. U' U = X, not U U'
  chol_cov <- chol(cov_mat)

  n_cols <- nrow(cov_mat)

  base_normal_mat <- matrix(rnorm(n_rows * n_cols),
                            nrow = n_rows,
                            ncol = n_cols
  )

  base_normal_mat %*% chol_cov
}
