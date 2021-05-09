#' predict_prepare
#'
#' Prepares design matrices and other inputs needed for predictive methods
#'
#' @param x_mat Original data, i.e. row for each unit.
#' @param index_list List of indices to use in the matching:
#'   \describe{
#'     \item{\code{treat_train}}{Index of treated units for the training matrix}
#'     \item{\code{control_train}}{Index of control units for the training matrix}
#'     \item{\code{treat_test}}{Index of treated units for the test matrix}
#'     \item{\code{control_test}}{Index of control units for the test matrix}
#'   }
#' @param design The design matrix form  to use.
#' @return list:
#'   \describe{
#'     \item{\code{x_train}}{Design matrix for training}
#'     \item{\code{x_test}}{Design matrix for testing}
#'     \item{\code{y_train}}{Training outcome vector}
#'     \item{\code{y_test}}{Test outcome vector}
#' }
#'
#' @export
predict_prepare <- function(x_mat,
                            index_list,
                            design = c(
                              "cross_all",
                              "cross_random",
                              "differences_random",
                              "differences_plain"
                            )) {
  design <- match.arg(design)

  stopifnot(length(index_list[["treat_train"]]) ==
              length(index_list[["control_train"]]) &&
              length(index_list[["treat_test"]]) ==
              length(index_list[["control_test"]]))

  ## ------------------------------------

  treat_train_ind <- index_list[["treat_train"]]
  control_train_ind <- index_list[["control_train"]]

  treat_test_ind <- index_list[["treat_test"]]
  control_test_ind <- index_list[["control_test"]]

  ## ------------------------------------

  x_treat_train <- x_mat[treat_train_ind, , drop = FALSE]
  x_control_train <- x_mat[control_train_ind, , drop = FALSE]

  x_treat_test <- x_mat[treat_test_ind, , drop = FALSE]
  x_control_test <- x_mat[control_test_ind, , drop = FALSE]

  train_index <- seq_len(length(treat_train_ind))
  test_index <- seq_len(length(treat_test_ind))

  ## y_train_full <-

  ## ------------------------------------

  if (design %in% c("cross_random", "differences_random")) {
    ## basically exactly half in each

    treat_left_train <- train_index %in%
      sample(length(treat_train_ind), size = length(treat_train_ind) / 2L)

    treat_left_test <- test_index %in%
      sample(length(treat_test_ind), size = length(treat_test_ind) / 2L)
  } else {
    ## cross_all has all, differences_random just has left
    treat_left_train <- rep(TRUE, times = length(treat_train_ind))
    treat_left_test <- rep(TRUE, times = length(treat_test_ind))
  }

  if (design == "cross_all") {
    treat_right_train <- rep(TRUE, times = length(treat_train_ind))
    treat_right_test <- rep(FALSE, times = length(treat_test_ind))

    train_order <- seq_len(length(treat_train_ind) * 2L)
    test_order <- seq_len(length(treat_test_ind))

    y_train <- rep(c(1L, 0L),
                   each = length(treat_train_ind)
    )
    y_test <- rep(1L, times = length(treat_test_ind))
  } else {
    treat_right_train <- !treat_left_train
    treat_right_test <- !treat_left_test

    train_order <- order(c(
      train_index[treat_left_train],
      train_index[treat_right_train]
    ))
    test_order <- order(c(
      test_index[treat_left_test],
      test_index[treat_right_test]
    ))

    y_train <- treat_left_train * 1L
    y_test <- treat_left_test * 1L
  }

  if (design %in% c("cross_all", "cross_random")) {
    ## these will "cross" the rows, and use all
    gen_design_mat <- function(x, y) {
      cbind(x, y)
    }
  } else {
    ## differences regime
    gen_design_mat <- function(x, y) {
      x - y
    }
  }

  ## ------------------------------------

  x_tc_train <- gen_design_mat(
    x_treat_train[treat_left_train, , drop = FALSE],
    x_control_train[treat_left_train, , drop = FALSE]
  )
  x_ct_train <- gen_design_mat(
    x_control_train[treat_right_train, , drop = FALSE],
    x_treat_train[treat_right_train, , drop = FALSE]
  )

  x_tc_test <- gen_design_mat(
    x_treat_test[treat_left_test, , drop = FALSE],
    x_control_test[treat_left_test, , drop = FALSE]
  )
  x_ct_test <- gen_design_mat(
    x_control_test[treat_right_test, , drop = FALSE],
    x_treat_test[treat_right_test, , drop = FALSE]
  )

  ## combining to form full input
  list(
    x_train = rbind(x_tc_train, x_ct_train)[train_order, , drop = FALSE],
    x_test = rbind(x_tc_test, x_ct_test)[test_order, , drop = FALSE],
    y_train = y_train,
    y_test = y_test
  )
}
