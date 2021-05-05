#' Run \code{compute_sim_result} in parallel using \code{parallel::mclapply}
#'
#' @inheritParams compute_sim_result
#' @param num_cores How many cores to use.
#' @param iterations How many iterations to do.
#' @param names_list ***
#'
#' @export
parallel_sim <- function(x_generator = default_x_generator,
                         treat_prob_generator,
                         mean_generator,
                         error_generator = default_error_generator,
                         n_sink_gen = n_sink_generator(),
                         match_method = "with_replacement",
                         n_rows = 500L,
                         n_cols = 5L,
                         num_weight_vectors = 100L,
                         num_cores = parallel::detectCores() - 1,
                         iterations = 100L,
                         names_list = NULL,
                         silent = !interactive()) {
  sims <- parallel::mclapply(1L:iterations, function(j) {
    if (!silent) {
      print(paste0("iteration ", j, "/", iterations))
    }
    compute_sim_result(
      x_generator = x_generator,
      treat_prob_generator = treat_prob_generator,
      mean_generator = mean_generator,
      error_generator = error_generator,
      n_sink_gen = n_sink_gen,
      match_method = match_method,
      n_rows = n_rows,
      n_cols = n_cols,
      num_weight_vectors = num_weight_vectors,
      silent = silent
    )
  }, mc.cores = num_cores)
}
