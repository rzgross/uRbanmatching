#' generate_plot_list
#'
#' Takes the result of \code{generate_p_cut_frame} and creates plot lists based on our choice of RMSE function.
#'
#' @param full_plot_frame Result from \code{generate_p_cut_frame}.
#' @param rmse_func Function that takes in a vector of estimates and computes an RMSE.
#' @param silent Whether to suppress messages as it's running.. Default \code{!interactive()}.
#'
#' @export
generate_plot_list <- function(full_plot_frame,
                               rmse_func = rmse_from_one_func,
                               silent = !interactive()) {
  rmse_boot <- function(vec) {
    sd(boot_func(vec, rmse_func, 500L))
  }

  n_row_split <- split(full_plot_frame, full_plot_frame$n_rows)

  lapply(n_row_split, function(df) {
    if (!silent) {
      print(df$n_rows[1])
    }
    df_by_p_cut <- split(df, df$p_cut)
    do.call(rbind, lapply(df_by_p_cut, function(x) {
      data.frame(
        p_cut = x$p_cut[1],
        rmse_naive = rmse_func(x$naive_est),
        se_naive = rmse_boot(x$naive_est),
        rmse_propensity = rmse_func(x$propensity_est),
        se_propensity = rmse_boot(x$propensity_est),
        rmse_mahal = rmse_func(x$mahal_est),
        se_mahal = rmse_boot(x$mahal_est),
        rmse_weighted = rmse_func(x$weighted_est),
        se_weighted = rmse_boot(x$weighted_est)
      )
    }))
  })
}
