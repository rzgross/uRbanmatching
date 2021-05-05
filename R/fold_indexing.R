#' Produces a list of indices split into k-folds
#'
#' @param index_length length of vector to index
#' @param num_folds how many folds you need
#' @return list of vectors of indices
#'
#' @export
fold_indexing <- function(index_length,
                          num_folds) {
  random_index <- sample(1L:index_length, size = index_length)

  gap_len <- index_length %% num_folds
  sub_length <- index_length - gap_len

  cut_indices <- seq(1L, sub_length + 1L,
                     length.out = num_folds + 1L
  ) +
    cumsum(c(0L, rep(1L, gap_len), rep(0L, num_folds - gap_len)))

  lapply(1L:num_folds, function(j) {
    random_index[cut_indices[j]:(cut_indices[j + 1L] - 1L)]
  })
}
