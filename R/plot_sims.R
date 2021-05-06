##' Plotting simulation results. Plan: four plots, one for each N value
##'
##' @param all_results result from \code{generate_plot_list}
##' @param y_rel where to limit y to above and below (length two vector)
##' @return Nothing, just plots to the device
##' @import grDevices
##' @import graphics
##' @export
plot_sims <- function(all_results,
                      y_rel = c(0.1, 0.4)) {
  x_push <- 0.15
  y_push <- 0.15

  plot(
    x = 0, y = 0,
    xlim = c(-x_push / 3, 2 + x_push * 4),
    ylim = c(0, 2 + y_push * 1),
    col = rgb(0, 0, 0, 0),
    xlab = "",
    ylab = "",
    axes = FALSE
  )

  par(xpd = TRUE)
  par(mar = c(1, 1, 1, 1))

  plot_results_block(all_results[["500"]],
                     x_lim = c(0, 0.95),
                     y_lim = c(.3, 1.8) + y_push,
                     y_rel
  )
  plot_results_block(all_results[["1000"]],
                     x_lim = c(1.05, 2) + x_push,
                     y_lim = c(.3, 1.8) + y_push,
                     y_rel
  )

  ## add legends, overall title
  text(
    x = c(0.475, 1.525) + x_push * c(0, 1),
    y = c(1.8, 1.8) + c(y_push, y_push),
    adj = c(0.5, 0.5),
    labels = paste("n = ", c(500, 1000)),
    cex = 1.2
  )

  y_seq <- seq(1.1 + y_push, 1.5 + y_push, length.out = 3)
  col_vec <- c(
    "red",
    "green",
    "purple"
  )

  segments(
    x0 = rep(2 + x_push * 2, 3),
    x1 = rep(2 + x_push * 2.8, 3),
    y0 = y_seq,
    col = col_vec,
    lty = c(1, 1, 1)
  )
  points(
    x = rep(2 + x_push * 2.4, 3),
    y = y_seq,
    pch = c(20, 20, 20),
    col = col_vec
  )
  text(
    x = rep(2 + x_push * 2.9, 3),
    y = y_seq,
    adj = c(0, 0.5),
    labels = c(
      "Propensity", "Mahalanobis", "Our Method"
    )
  )

  rect(
    xleft = 2 + x_push * 2.385,
    xright = 2 + x_push * 2.415,
    ybottom = 0.8,
    ytop = 1.05,
    col = 1, border = NA
  )
  points(
    x = 2 + x_push * 2.4,
    y = 0.925,
    pch = 20,
    cex = 1.8
  )

  text(
    x = rep(2 + x_push * 2.9, 2),
    y = c(0.95, 0.88),
    adj = c(0, 0.5),
    cex = c(1, 0.8),
    labels = c("CI, +-2 se", "(bootstrapped)")
  )
}
