#' brier_score
#'
#' Takes in training/test data and a prediction function to use, then generates brier score.
#'
#' @param train_test_list Output from \code{predict_prepare}
#' @param match_predict_function Function to predict treated unit, either  \code{match_predict_xgb} or \code{match_predict_linear}
#' @param avg Whether to average (or sum) the briers.
#' @return Result from \code{calc_brier}, length-one double
#'
#' @export
brier_score <- function(train_test_list,
                        match_predict_function = match_predict_xgb(),
                        avg = TRUE) {
  calc_brier(match_predict_function(train_test_list),
             train_test_list[["y_test"]],
             avg = avg
  )
}
