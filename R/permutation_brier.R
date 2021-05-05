#' For a given match, computes the brier score distribution
#' if the pairing were truly random
#'
#' Works for bipartite and non-bipartite, since it just
#' switches the labels within pairs.
#' @param x_mat typical input matrix
#' @param match_list match result
#' @param design see \code{predict_prepare}
#' @param use_cv logical, default TRUE: use CV to get briers? Else
#'   split.
#' @param num_permutations how many permutations to do
#' @param match_predict_function function to predict treated units
#' @param num_folds if using CV, how many folds?
#' @param train_fraction if using split, fraction to train?
#' @return vector of brier scores for random pairings
#'
#' @export
permutation_brier <- function(x_mat,
                              match_list,
                              design = "cross_all",
                              use_cv = TRUE,
                              num_permutations = 100L,
                              match_predict_function = match_predict_xgb(),
                              num_folds = 5,
                              train_fraction = 0.7) {
  if (use_cv) {
    if (!missing(train_fraction)) {
      stop(
        "only set `train_fraction` if not using cross-validation ",
        "(set `use_cv = FALSE`)"
      )
    }

    brier_function <- (function(x_mat,
                                design,
                                num_folds,
                                match_predict_function) {
      function(match_list) {
        brier_score_cv(
          x_mat,
          match_list,
          design,
          num_folds,
          match_predict_function
        )
      }
    })(x_mat, design, num_folds, match_predict_function)
  } else {
    if (!missing(num_folds)) {
      stop(
        "only set `num_folds` if using cross-validation ",
        "(set `use_cv = TRUE`)"
      )
    }

    brier_function <- (function(x_mat,
                                design,
                                num_folds,
                                match_predict_function) {
      function(match_list) {
        brier_score_split(
          x_mat,
          match_list,
          design,
          train_fraction,
          match_predict_function
        )
      }
    })(x_mat, design, train_fraction, match_predict_function)
  }

  ## ------------------------------------

  full_index <- seq_len(length(match_list[["treat_index"]]))

  unlist(lapply(1L:num_permutations, function(i) {
    swaps <- full_index %in% sample(full_index,
                                    size = floor(length(full_index) / 2)
    )

    swap_match <- swap_pairs(
      match_list,
      swaps
    )

    brier_function(swap_match)
  }))
}
