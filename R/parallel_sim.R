#' parallel_sim
#'
#' Run \code{compute_sim_result} in parallel using \code{parallel::mclapply} or \code{parallel::parLapply}
#'
#' @inheritParams compute_sim_result
#' @param num_cores Number of cores to use.
#' @param num_cores1 Number of cores to use.
#' @param iterations Number of iterations to complete.
#' @param names_list Default NULL. List of index names to return with the results.
#' @param operating_system "Mac," "Windows," etc. If Windows, uses parLapply instead of mclapply.
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
                         iterations = 100L,
                         num_cores = parallel::detectCores(),
                         num_cores1 = 1,
                         operating_system = "Mac",
                         names_list = NULL,
                         silent = !interactive()) {
  if(operating_system != "Windows"){
    sims <- parallel::mclapply(1L:max(1,iterations-1), function(j) {
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
    }, mc.cores = num_cores1)
  }
  else{
    sims <- parallel::parLapply(1L:max(1,iterations-2), function(j) {
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
    }, mc.cores = num_cores1)
  }
}

