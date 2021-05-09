#' match_predict_xgb
#'
#' Function to predict treatment / control pairs using xgboost. The returned function takes in training and test data (output from \code{predict_prepare}), trains an xgboost model on the training, predicts on the test, and returns the test vector
#'
#' @param nrounds Training rounds for the xgb algorithm.
#' @param nthread Number of threads to use for fitting, default 1.
#' @param params List of params to pass to xgboost, most likely something like \code{eta} and \code{max.depth}
#' @return Returns a function that takes in a \code{train_test_list} from \code{predict_prepare}; this function returns a vector of predictions for the test data.
#' @param ... Additional xgboost params.
#' @export
match_predict_xgb <- function(nrounds = 50,
                              nthread = 1,
                              params = list(eta = 0.1, max.depth = 4),
                              ...) {
  if (!requireNamespace("xgboost")) {
    stop("xgboost not installed", call. = FALSE)
  }

  function(train_test_list) {
    xgb_train <- xgboost::xgb.DMatrix(
      train_test_list[["x_train"]],
      label = train_test_list[["y_train"]]
    )

    xgb_test <- xgboost::xgb.DMatrix(
      train_test_list[["x_test"]],
      label = train_test_list[["y_test"]]
    )

    train_res <- xgboost::xgb.train(
      params = params,
      data = xgb_train,
      nrounds = nrounds,
      nthread = nthread,
      objective = "binary:logistic",
      ...
    )

    predict(train_res, newdata = xgb_test)
  }
}
