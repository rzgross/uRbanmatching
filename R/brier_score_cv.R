#' For a match, calculates brier score using cross validation
#'
#' @inheritParams permutation_brier
#' @param num_folds how many CV folds to use
#' @return brier score (averaged over all units)
#'
#' @export
brier_score_cv <- function(x_mat,
                           match_list,
                           design = "cross_all",
                           num_folds = 5,
                           match_predict_function = match_predict_xgb()) {
  k_fold_lists <- generate_k_fold_index(
    match_list,
    num_folds
  )

  pred_list <- lapply(k_fold_lists, function(k_fold) {
    train_test_list <- predict_prepare(x_mat,
                                       k_fold,
                                       design = design
    )
    c(
      brier_score(train_test_list,
                  match_predict_function,
                  avg = FALSE
      ),
      length(train_test_list[["y_test"]])
    )
  })

  sums <- Reduce("+", pred_list)
  sums[1] / sums[2]
}
