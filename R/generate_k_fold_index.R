#' Constructs a k-fold list of \code{index_list} objects for a given
#' \code{match_list}
#'
#' @param match_list typical \code{match_list} entry
#' @param num_folds how many folds you want, default 5
#'
#' @export
generate_k_fold_index <- function(match_list,
                                  num_folds = 5L) {
  length_index <- length(match_list[["treat_index"]])
  all_index <- 1L:length_index

  fold_res <- fold_indexing(length_index, num_folds)

  lapply(fold_res, function(inds) {
    index_list_from_match(
      match_list,
      !(all_index %in% inds)
    )
  })
}
