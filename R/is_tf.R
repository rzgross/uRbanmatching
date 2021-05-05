#' Test if \code{x} is a length-one logical
#'
#' @param x ideally a logical!
#' @return TRUE if x is a length one logical, else FALSE
#'
#' @keywords internal
is_tf <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}
