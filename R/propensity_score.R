#' propensity_score
#'
#' Calculates propensity scores for a given matrix and treatment vector. Takes in an input matrix and a treatment vector, along with a function that makes predictions (default xgboost method given) and returns a predicted probability of treatment for each unit, either using in-sample or out-of-sample fits.
#' @param x_mat Standard input matrix (already rank adjusted).
#' @param treat_vec Usual 0/1 treatment vector.
#' @param propensity_list See \code{gen_propensity_list}
#' @return Returns a vector equal in length to treat_vec of propensity score.
#'
#' @export
propensity_score <- function(x_mat,
                             treat_vec,
                             propensity_list = gen_propensity_list()) {
  propensity_function <- propensity_list[["propensity_function"]]

  if (!propensity_list[["oos_propensity"]]) {
    ## very simple...:
    train_treat_list <- list(
      x_train = x_mat,
      x_test = x_mat,
      y_train = treat_vec,
      y_test = treat_vec
    )
    return(propensity_function(train_treat_list))
  }

  ## ------------------------------------

  fold_res <- fold_indexing(nrow(x_mat), propensity_list[["n_folds"]])

  prediction_result <- lapply(fold_res, function(fold_inds) {
    fold_lgl <- (1L:nrow(x_mat)) %in% fold_inds
    train_treat_list <- list(
      x_train = x_mat[!fold_lgl, , drop = FALSE],
      x_test = x_mat[fold_inds, , drop = FALSE],
      y_train = treat_vec[!fold_lgl],
      y_test = treat_vec[fold_inds]
    )
    propensity_function(train_treat_list)
  })

  unlist(prediction_result)[order(unlist(fold_res))]
}
