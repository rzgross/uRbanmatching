#' generate_p_cut_frame
#'
#' Processes the data into the format required for plotting.
#'
#' @param full_raw_results Dataframe of results (spec to come)
#' @param p_cut_vals Vector of cut values to plot.
#' @return List per \code{n_rows} value of a dataframe ready to be plotted.
#' @param silent Whether to suppress messages as it's running.. Default \code{!interactive()}.
#'
#' @export
generate_p_cut_frame <- function(full_raw_results,
                                 p_cut_vals = seq(from = 0.2, to = 1, by = 0.05),
                                 silent = !interactive()) {
  per_model <- split(full_raw_results, full_raw_results$id)

  p_cut_sims <- lapply(p_cut_vals, function(p_cut) {
    if (!silent) {
      print(p_cut)
    }
    sims_frame <- do.call(rbind, lapply(per_model, function(df) {
      df_filter <- df[df$p_brier <= p_cut, ]
      if (nrow(df_filter) > 0) {
        ## we increase sinks till we're done,
        ## so take the one with the least sinks
        return(df_filter[which.min(df_filter$n_sinks), ])
      }
      if (all(df$p_brier > 0.99999)) { # safely '1'-equality checking
        ## This is very likely the same as argmax sinks
        return(df[which.max(df$raw_brier), ])
      }
      return(df[which.min(df$p_brier), ])
    }))
    sims_frame$too_big <- sims_frame$p_brier > p_cut
    sims_frame$p_cut <- p_cut
    sims_frame
  })

  do.call(rbind, p_cut_sims)
}
