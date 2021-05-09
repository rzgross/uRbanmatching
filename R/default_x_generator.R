#' default_x_generator
#'
#' Function that generates an input data matrix. Used as default in \code{generate_simulation_input}.
#'
#' @param n_rows Number ofrows to produce.
#' @param n_cols Number of columns to produce.
#' @return A matrix of data.
#'
#' @export
default_x_generator <- function(n_rows,
                                n_cols) {
  sig_mat_pre <- matrix(
    stats::runif((n_cols - 1) * n_cols, 0, 0.1),
    n_cols, n_cols - 1
  )
  sig_mat_cov <- sig_mat_pre %*% t(sig_mat_pre) +
    diag(x = stats::runif(n_cols, 0, 0.3))
  sig_mat_cor <- stats::cov2cor(sig_mat_cov)

  x_from_cov(sig_mat_cor,
             n_rows = n_rows
  )
}
