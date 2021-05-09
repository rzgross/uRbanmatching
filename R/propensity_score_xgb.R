#' propensity_score_xgb
#'
#' Function factory to predict treatment using xgboost.
#'
#' @inheritParams match_predict_xgb
#' @return returns a function that accepts \code{train_test_list} and this returns a vector of predictions for the test data
#'
#' @export
propensity_score_xgb <- function(nrounds = 50,
                                 nthread = 1,
                                 params = list(eta = 0.1, max.depth = 3),
                                 ...) {
  match_predict_xgb(nrounds = nrounds,
                    nthread = nthread,
                    params = params,
                    ...)
}
