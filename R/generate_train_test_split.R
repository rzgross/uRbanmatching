#' Creates an \code{index_list} from a \code{match_list}, splitting
#' according to \code{train_fraction}
#'
#' @param match_list typical \code{match_list} entry
#' @param train_fraction fraction (between 0 and 1) to
#'   use for training data (and the rest for test)
#'
#' @export
generate_train_test_split <- function(match_list,
                                      train_fraction = 0.7) {
  stopifnot(train_fraction >= 0 && train_fraction <= 1)

  length_index <- length(match_list[["treat_index"]])
  all_index <- 1L:length_index

  train_index <- all_index %in%
    sample(all_index, size = floor(train_fraction * length_index))
  index_list_from_match(
    match_list,
    train_index
  )
}
