#' For a match, calculates brier score on a test split of the data
#'
#' @inheritParams permutation_brier
#' @param train_fraction split of data to use for training
#' @return brier score
#'
#' @export
brier_score_split <- function(x_mat,
                              match_list,
                              design = "cross_all",
                              train_fraction = 0.7,
                              match_predict_function = match_predict_xgb()) {
  brier_score(predict_prepare(
    x_mat,
    generate_train_test_split(match_list, train_fraction),
    design
  ),
  match_predict_function,
  avg = TRUE
  )
}
