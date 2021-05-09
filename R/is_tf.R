#' is_tf
#'
#' Test if \code{x} is a length-one logical.
#'
#' @param x Object to test.
#' @return TRUE if x is a length one logical, else FALSE
#'
#' @export
is_tf <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}
