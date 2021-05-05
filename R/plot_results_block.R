plot_results_block <- function(n_results,
                               x_lim,
                               y_lim,
                               y_rel = c(0, 1),
                               rect_width = 0.006,
                               x_shift = 0.015) {
  n_results$high_naive <- n_results$rmse_naive +
    2 * n_results$se_naive
  n_results$low_naive <- pmax(n_results$rmse_naive -
                                2 * n_results$se_naive, 0)
  n_results$high_propensity <- n_results$rmse_propensity +
    2 * n_results$se_propensity
  n_results$low_propensity <- pmax(n_results$rmse_propensity -
                                     2 * n_results$se_propensity, 0)
  n_results$high_mahal <- n_results$rmse_mahal + 2 * n_results$se_mahal
  n_results$low_mahal <- pmax(n_results$rmse_mahal - 2 * n_results$se_mahal, 0)
  n_results$high_weighted <- n_results$rmse_weighted + 2 * n_results$se_weighted
  n_results$low_weighted <- pmax(n_results$rmse_weighted - 2 * n_results$se_weighted, 0)

  max_height <- y_rel[2]

  ## max_height <- min(
  ##     max(
  ##         max(n_results[["high_propensity"]]),
  ##         max(n_results[["high_mahal"]]),
  ##         max(n_results[["high_weighted"]])
  ##     ) * 1.05,
  ##     y_rel
  ## )
  min_x <- min(n_results[["p_cut"]])

  x_adj <- function(x) {
    ((x - min_x) / (1 - min_x)) * (x_lim[2] - x_lim[1]) + x_lim[1]
  }
  y_adj <- function(y) {
    (y / max_height) * (y_lim[2] - y_lim[1]) + y_lim[1]
  }

  rect_width <- rect_width * (1 - min_x)
  x_shift <- x_shift * (1 - min_x)
  x_ax_adj <- 0.05 * (1 - min_x)
  y_ax_adj <- 0.02 * max_height

  ## ------------------------------------

  ## axes
  segments(
    x0 = x_adj(min_x - x_ax_adj),
    x1 = x_adj(1),
    y0 = y_adj(0)
  )
  x_tick_pre_seq <- unique(n_results[["p_cut"]])
  x_tick_seq <- x_adj(x_tick_pre_seq)
  segments(
    x0 = x_tick_seq,
    y0 = y_adj(-y_ax_adj),
    y1 = y_adj(0)
  )
  text(
    x = x_tick_seq,
    y = y_adj(-y_ax_adj * 2),
    labels = round(x_tick_pre_seq, 2),
    adj = c(0.5, 0.5),
    cex = 0.5
  )

  text(
    x = x_adj((min_x + 1) / 2),
    y = y_adj(-y_ax_adj * 5),
    labels = expression("p"["cut"]),
    adj = c(0.5, 0.5),
    cex = 1
  )

  segments(
    x0 = x_adj(min_x - x_ax_adj),
    y0 = y_adj(0),
    y1 = y_adj(max_height)
  )
  y_tick_pre_seq <- seq(
    from = 0,
    to = max_height,
    length.out = 5
  )
  y_tick_seq <- y_adj(y_tick_pre_seq)
  segments(
    y0 = y_tick_seq,
    x0 = x_adj(min_x - x_ax_adj * 1.5),
    x1 = x_adj(min_x - x_ax_adj)
  )
  text(
    x = x_adj(min_x - x_ax_adj * 2),
    y = y_tick_seq,
    labels = round(y_tick_pre_seq, 2),
    adj = c(0.5, 0.5),
    cex = 0.5
  )

  text(
    x = x_adj(min_x - x_ax_adj * 3),
    y = y_adj(max_height / 2),
    labels = "RMSE",
    adj = c(0.5, 0.5),
    cex = 1,
    srt = 90
  )

  par(xpd = FALSE)

  ## ------------------------------------
  ## naive, may not even show!

  ## naive_mean <- n_results[["rmse_naive"]][1]
  ## naive_high <- n_results[["high_naive"]][1]
  ## naive_low <- n_results[["low_naive"]][1]

  ## if (naive_mean / max_height < 1.1) {
  ##     segments(
  ##         x0 = x_lim[1],
  ##         x1 = x_lim[2],
  ##         y0 = y_adj(naive_mean),
  ##         lty = 2, col = 1, lwd = 1
  ##     )
  ##     segments(
  ##         x0 = x_lim[1],
  ##         x1 = x_lim[2],
  ##         y0 = y_adj(naive_high),
  ##         lwd = 0.7,
  ##         lty = 2, col = rgb(0, 0, 0, 0.5)
  ##     )
  ##     segments(
  ##         x0 = x_lim[1],
  ##         x1 = x_lim[2],
  ##         y0 = y_adj(naive_low),
  ##         lwd = 0.7,
  ##         lty = 2, col = rgb(0, 0, 0, 0.5)
  ##     )
  ## }

  ## ------------------------------------
  ## propensity

  points(x_adj(n_results[["p_cut"]] - x_shift),
         y_adj(n_results[["rmse_propensity"]]),
         type = "l", col = rgb(0.7, 0.3, 0.2, 1)
  )
  points(x_adj(n_results[["p_cut"]] - x_shift),
         y_adj(n_results[["rmse_propensity"]]),
         pch = 20,
         col = rgb(0.7, 0.3, 0.2, 1)
  )
  rect(
    xleft = x_adj(n_results[["p_cut"]] - rect_width - x_shift),
    xright = x_adj(n_results[["p_cut"]] + rect_width - x_shift),
    ybottom = y_adj(n_results[["low_propensity"]]),
    ytop = y_adj(n_results[["high_propensity"]]),
    col = rgb(0.7, 0.3, 0.2, 0.8),
    border = NA
  )

  ## ------------------------------------
  ## mahal

  points(x_adj(n_results[["p_cut"]] + x_shift),
         y_adj(n_results[["rmse_mahal"]]),
         type = "l", col = rgb(0.2, 0.8, 0.4, 1)
  )
  points(x_adj(n_results[["p_cut"]] + x_shift),
         y_adj(n_results[["rmse_mahal"]]),
         pch = 20,
         col = rgb(0.2, 0.8, 0.4, 1)
  )
  rect(
    xleft = x_adj(n_results[["p_cut"]] - rect_width + x_shift),
    xright = x_adj(n_results[["p_cut"]] + rect_width + x_shift),
    ybottom = y_adj(n_results[["low_mahal"]]),
    ytop = y_adj(n_results[["high_mahal"]]),
    col = rgb(0.2, 0.8, 0.4, 0.8),
    border = NA
  )

  ## ------------------------------------
  ## weighted

  points(x_adj(n_results[["p_cut"]]),
         y_adj(n_results[["rmse_weighted"]]),
         type = "l", col = rgb(0.3, 0.1, 0.9, 1)
  )
  points(x_adj(n_results[["p_cut"]]),
         y_adj(n_results[["rmse_weighted"]]),
         pch = 20,
         col = rgb(0.3, 0.1, 0.9, 1)
  )
  rect(
    xleft = x_adj(n_results[["p_cut"]] - rect_width),
    xright = x_adj(n_results[["p_cut"]] + rect_width),
    ybottom = y_adj(n_results[["low_weighted"]]),
    ytop = y_adj(n_results[["high_weighted"]]),
    col = rgb(0.3, 0.1, 0.9, 0.8),
    border = NA
  )

  par(xpd = TRUE)

  rect(
    xleft = x_adj(min_x - 3 * rect_width - x_shift),
    xright = x_adj(1.1),
    ybottom = y_adj(y_rel[2] * 1.05),
    ytop = 50,
    border = NA, col = rgb(1, 1, 1, 1)
  )

}
