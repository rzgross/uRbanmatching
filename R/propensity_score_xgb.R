#' Function factory to predict treatment using xgboost
#'
#' This operates the exact same as \code{match_predict_xgb}, which
#' is itself quite a simple wrap of regular xgboost code.
#' Main difference here is the input data is even simpler.
#' The returned function accepts one parameter, \code{train_test_list},
#' a list with \code{x_train}, \code{x_test}, \code{y_train}, \code{y_test}
#' @inheritParams match_predict_xgb
#' @return returns a function that accepts \code{train_test_list}
#' and this returns a vector of predictions for the test data
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
