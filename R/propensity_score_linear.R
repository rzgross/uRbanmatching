#' Function factory to predict treatment using \code{glm} (binomial)
#' or \code{lm}
#'
#' Does simple wrap around \code{glm} / \code{lm}
#' The returned function accepts one parameter, \code{train_test_list},
#' a list with \code{x_train}, \code{x_test}, \code{y_train}, \code{y_test}
#' @return returns a function that accepts \code{train_test_list}
#' and this returns a vector of predictions for the test data
#' @param use_linear_lm ***
#' @export
propensity_score_linear <- function(use_linear_lm = FALSE) {
  function(train_test_list) {
    train_frame <- as.data.frame(train_test_list[["x_train"]])
    train_frame[["y"]] <- train_test_list[["y_train"]]

    test_frame <- as.data.frame(train_test_list[["x_test"]])

    if (use_linear_lm) {
      train_res <- lm(y ~ ., data = train_frame)
      lin_pred <- predict(train_res, newdata = test_frame,
                          type = "response")
      return(pmax(pmin(lin_pred, 1), 0))
    }

    train_res <- glm(y ~ ., data = train_frame, family = "binomial")
    predict(train_res, newdata = test_frame, type = "response")
  }
}
