#' permutation_brier
#'
#' For a given match, computes the brier score distribution, assuming the pairings were truly random.
#'
#' @param x_mat Typical input matrix
#' @param match_list Match result
#' @param design See \code{predict_prepare}
#' @param use_cv Logical, default TRUE: wether to use CV to get briers? Else split.
#' @param num_permutations Number of permutations to do.
#' @param match_predict_function Function to predict treated units.
#' @param num_folds If using CV, how many folds to make.
#' @param train_fraction If using split, what fraction to train.
#' @return Vector of brier scores for random pairings
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
