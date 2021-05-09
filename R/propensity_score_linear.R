#' propensity_score_linear
#'
#' Function to predict treatment using \code{glm} (binomial) or \code{lm}.
#'
#' @param use_linear_lm Whether to use lm or glm.
#'
#' @return returns a function that accepts \code{train_test_list}
#' and this returns a vector of predictions for the test data
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
