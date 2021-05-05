#' Takes in training/test data and a prediction
#' function to use, generates brier score
#'
#' @param train_test_list output from \code{predict_prepare}
#' @param match_predict_function function to predict treated unit,
#'   see e.g. \code{match_predict_xgb}
#'   (produces the default on calling) or
#'   \code{match_predict_linear}
#' @param avg logical, default TRUE: should we average (or sum) the briers?
#' @return result from \code{calc_brier}, length-one double
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
