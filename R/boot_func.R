#' boot_func
#'
#' Resampling bootstrap of a vector for a given function.
#'
#' @param vec Vector to bootstrap resample from.
#' @param func Function to call on the samples.
#' @param resamples How many resamples to do.
#' @return Returns vector of length \code{resamples}, each from \code{func} applied to a bootstrap resample.
#' @export

boot_func <- function(vec, func, resamples = 100L) {
  rep_list <- lapply(1L:resamples, function(j) {
    func(sample(vec, replace = TRUE))
  })
  unlist(rep_list)
}
