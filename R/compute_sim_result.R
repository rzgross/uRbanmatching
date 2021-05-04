#' Takes in functions to generate simulation data, and computes
#' simulation results for our method and some competitors
#'
#' @inheritParams generate_simulation_input
#' @inheritParams all_bipartite_matches
#' @param n_sink_gen Default \code{n_sink_generator}, and you'll
#'   probably want to use that: this argument should be a function that
#'   accepts a \code{treat_vec} and produces a vector of sink numbers.
#' @param num_weight_vectors How many weight vectors to generate.
#' @param silent Default \code{!interactive()}, if you want to
#'   suppress messages.
#' @return Returns a named list:
#' \describe{
#'     \item{\code{naive_est}}{Just a number: mean difference between all
#'   treated units and all control}
#'     \item{\code{propensity_results}}{List of lists: each with
#'   \code{n_sinks} and the \code{est}}
#'     \item{\code{mahal_results}}{Same as above}
#'     \item{\code{weighted_results}}{List of lists: each with
#'   \code{n_sinks}, the raw brier score, the permutation brier score,
#'   and the \code{est}}
#' }
#'
#' @export
compute_sim_result <- function(x_generator = default_x_generator,
                               treat_prob_generator,
                               mean_generator,
                               error_generator = default_error_generator,
                               n_sink_gen = n_sink_generator(),
                               match_method = "with_replacement",
                               n_rows = 500L,
                               n_cols = 5L,
                               num_weight_vectors = 100L,
                               silent = !interactive()) {
  sim_data <- generate_simulation_input(
    n_rows = n_rows,
    n_cols = n_cols,
    x_generator = x_generator,
    treat_prob_generator = treat_prob_generator,
    mean_generator = mean_generator,
    error_generator = error_generator
  )

  x_mat <- sim_data[["x_mat"]]
  y_vector <- sim_data[["y_vec"]]
  treat_vec <- sim_data[["treat_vec"]]

  rm(sim_data)

  n_sinks <- n_sink_gen(treat_vec)
  match_list_est_func <- (function(y_vector, treat_vec) {
    function(match_list) {
      ## TODO: this should be regression eval
      ## leaving here to make current simulations replicable
      match_estimate(
        match_list = match_list,
        y_vector = y_vector,
        treat_vec = treat_vec
      )
    }
  })(y_vector, treat_vec)

  list_est_func <- (function(n_sinks) {
    function(match_lists) {
      Map(
        function(n_sink, match_list) {
          list(
            n_sinks = n_sink,
            est = match_list_est_func(match_list)
          )
        },
        n_sinks,
        match_lists
      )
    }
  })(n_sinks)

  naive_est <- mean(y_vector[treat_vec == 1]) - mean(y_vector[treat_vec == 0])

  ## ------------------------------------

  if (!silent) {
    message("propensity matches")
  }

  propensity_matches <- propensity_bipartite_matches(
    x_mat = x_mat,
    treat_vec = treat_vec,
    match_method = match_method,
    propensity_list = gen_propensity_list(
      propensity_function = propensity_score_linear(),
      oos_propensity = FALSE
    ),
    n_sinks = n_sinks
  )

  propensity_ests <- list_est_func(propensity_matches)

  ## ------------------------------------

  if (!silent) {
    message("mahal matches")
  }

  mahal_matches <- all_bipartite_matches(
    x_mat = x_mat,
    cov_x = covariance_with_ranks(x_mat),
    weight_list = list(rep(1 / n_cols, times = n_cols)),
    treat_vec = treat_vec,
    match_method = match_method,
    n_sinks = n_sinks
  )[[1]]

  mahal_ests <- list_est_func(mahal_matches)

  ## ------------------------------------

  if (!silent) {
    message("all weighted matches")
  }

  weight_list <- generate_random_weights(
    prior_weights = rep(1 / n_cols, times = n_cols),
    number_vectors = num_weight_vectors,
    minimum_weights = rep(1 / (3 * n_cols), times = n_cols)
  )

  brier_matches <- brier_bipartite_matches(
    x_mat = x_mat,
    cov_x = covariance_with_ranks(x_mat),
    weight_list = weight_list,
    treat_vec = treat_vec,
    match_method = match_method,
    n_sinks = n_sinks,
    silent = silent
  )

  permutation_results <- permutation_matches(
    matches_by_sinks = brier_matches[["matches_by_sinks"]],
    briers_by_sinks = brier_matches[["briers_by_sinks"]],
    x_mat = x_mat,
    n_sinks = n_sinks,
    silent = silent
  )

  weighted_results <- lapply(
    permutation_results[["best_matches"]],
    function(match_results) {
      list(
        n_sinks = match_results[["n_sinks"]],
        raw_brier = match_results[["raw_brier"]],
        permutation_brier = match_results[["permutation_brier"]],
        est = match_list_est_func(match_results[["match_list"]])
      )
    }
  )

  ## ------------------------------------

  list(
    naive_est = naive_est,
    propensity_results = propensity_ests,
    mahal_results = mahal_ests,
    weighted_results = weighted_results
  )
}


#' Takes in elements needed for propensity work, checks input and
#' builds a named list.
#'
#' @param propensity_function A function that accepts a list
#'   with four elements: \code{x_train}, \code{x_test},
#'   \code{y_train}, \code{y_test}, and forms predictions
#'   using \code{x_test} (I guess \code{y_test} isn't used yet)
#' @param oos_propensity Logical, do you want to predict out of sample
#'   for the propensity score? Most people don't, and indeed \code{FALSE}
#'   is the default.
#' @param n_folds Default NULL; how many folds you want if using
#'   out of sample propensity.
#' @return Named list, same names as input params.
#'
#' @export
gen_propensity_list <- function(propensity_function = propensity_score_xgb(),
                                oos_propensity = FALSE,
                                n_folds = NULL) {
  stopifnot(is_tf(oos_propensity))

  if (missing(oos_propensity) && !is.null(n_folds)) {
    oos_propensity <- TRUE
  }

  if (!oos_propensity) {
    if (!is.null(n_folds)) {
      stop("n_folds shouldn't be set if not using out of sample")
    }
  } else {
    n_folds <- ifelse(is.null(n_folds), 5L, n_folds)
    stopifnot(is.numeric(n_folds) && length(n_folds) == 1L &&
                !is.na(n_folds))
  }

  list(
    propensity_function = propensity_function,
    oos_propensity = oos_propensity,
    n_folds = n_folds
  )
}


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


#' Function factory to predict treatment using \code{glm} (binomial)
#' or \code{lm}
#'
#' Does simple wrap around \code{glm} / \code{lm}
#' The returned function accepts one parameter, \code{train_test_list},
#' a list with \code{x_train}, \code{x_test}, \code{y_train}, \code{y_test}
#' @inheritParams propensity_score_linear
#' @return returns a function that accepts \code{train_test_list}
#' and this returns a vector of predictions for the test data
#'
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


#' computes the covariance of the input matrix, adjusting for ranks if given
#'
#' adjusts to ranks if given; fixes ties; rescales etc;
#' returns the covariance of x_mat. If no rank_cols are given, this
#' function is identical to just calling cov(x_mat)
#' @param x_mat matrix of variables (numeric)
#' @param rank_cols names of columns to be converted to ranks before analysis
#' @return covariance matrix of x_mat (potentially adjusted for ranks)
#'
#' @export
covariance_with_ranks <- function(x_mat,
                                  rank_cols = NULL) {
  if (!is.matrix(x_mat)) {
    stop("x_mat must be a matrix")
  }

  ## convert to ranks
  x_mat <- ranked_x(
    x_mat,
    rank_cols
  )

  cov_x <- cov(x_mat)

  if (is.null(rank_cols)) {
    return(cov_x)
  }

  rank_cols <- rank_integer_index(rank_cols, x_mat)

  ## fix the rank ties variances
  var_untied <- var((1:nrow(x_mat)) / nrow(x_mat))
  sd_ratio <- sqrt(var_untied / diag(cov_x))

  ## don't fix the ones that aren't ranked
  sd_ratio[!(1L:ncol(x_mat) %in% rank_cols)] <- 1

  diag(sd_ratio) %*% cov_x %*% diag(sd_ratio)
}


#' Computes all matches, then gets the brier scores for each. Reorder by
#' number of sinks.
#'
#' @inheritParams all_bipartite_matches
#' @inheritParams brier_score_cv
#' @return List of matches within sink values,
#'  and brier scores for each.
#'
#' @export
brier_bipartite_matches <- function(x_mat,
                                    cov_x,
                                    weight_list,
                                    treat_vec,
                                    match_method = c(
                                      "with_replacement",
                                      "optimal",
                                      "greedy"
                                    ),
                                    n_sinks = 0L,
                                    caliper_list = gen_caliper_list(),
                                    propensity_list =
                                      match_propensity_list(NULL),
                                    sqrt_mahal = TRUE,
                                    tol_val = NULL,
                                    design = "cross_all",
                                    num_folds = 5,
                                    match_predict_function =
                                      match_predict_xgb(),
                                    silent = !interactive()) {
  if (is.null(n_sinks)) {
    n_sinks <- 0L
  }

  ## generate all matches: one per weight vector per n_sink value
  all_matches <- all_bipartite_matches(
    x_mat = x_mat,
    cov_x = cov_x,
    weight_list = weight_list,
    treat_vec = treat_vec,
    match_method = match_method,
    n_sinks = n_sinks,
    caliper_list = caliper_list,
    propensity_list = propensity_list,
    sqrt_mahal = sqrt_mahal,
    tol_val = tol_val
  )

  if (!silent) {
    message("getting briers")
  }

  ## get all brier scores for all results
  briers_by_sinks <- lapply(all_matches, function(all_by_sink) {
    if (!silent) {
      print(all_by_sink[[1]]["num_sinks"])
    }
    unlist(lapply(all_by_sink, function(indiv_match_list) {
      brier_score_cv(
        x_mat = x_mat,
        match_list = indiv_match_list,
        design = design,
        num_folds = num_folds,
        match_predict_function = match_predict_function
      )
    }))
  })

  list(
    matches_by_sinks = all_matches,
    briers_by_sinks = briers_by_sinks
  )
}



#' Takes matches and their Brier scores, and computes
#' permutation Brier scores and the best matches
#'
#' Can work for bipartite or non-bipartite matches,
#' the permutation is just over the labels.
#' @param matches_by_sinks List by number of sinks, each a list
#'   of match results (a match list), for each weight vector.
#' @param briers_by_sinks List by number of sinks, each a vector
#'   of Brier results. Basically a number for each match in
#'   \code{matches_by_sinks}.
#' @param x_mat Typical input matrix
#' @param n_sinks Vector of number of sinks - probably could get this
#'   directly from \code{matches_by_sinks}, but nice to be explicit.
#' @param approximate_by_best Logical, default \code{TRUE}. Only
#'   compute one permutation distribution, using the best result by brier
#'   score to do so. Useful because it changes little, but saves a ton of
#'   time.
#' @param silent Do you want to suppress message output? Default
#'   \code{!interactive()}.
#' @return Returns a list of two lists. The first is vectors
#'   of permutation Brier scores (one per match). The second is the
#'   best match at each sink value, along with some extra info about that
#'   match.
#' @author Colman Humphrey
#'
#' @export
permutation_matches <- function(matches_by_sinks,
                                briers_by_sinks,
                                x_mat,
                                n_sinks = 0L,
                                approximate_by_best = TRUE,
                                silent = !interactive()) {
  if (is.null(n_sinks)) {
    n_sinks <- 0L
  }

  if (length(n_sinks) != length(matches_by_sinks)) {
    stop("`n_sinks` is length ", length(n_sinks),
         " but `matches_by_sinks` is length ", length(matches_by_sinks))
  }

  if (length(n_sinks) != length(briers_by_sinks)) {
    stop("`n_sinks` is length ", length(n_sinks),
         " but `briers_by_sinks` is length ", length(briers_by_sinks))
  }

  best_brier_inds <- lapply(briers_by_sinks, function(x) {
    which(rank(x, ties.method = "first") == length(x))
  })

  if (!silent) {
    message("running permutations, will be a little slow")
  }

  permutation_briers <- lapply(seq_len(length(n_sinks)), function(j) {
    if (!silent) {
      message(paste0("Running permutations for ",
                     n_sinks[j], " sinks"))
    }
    if (approximate_by_best) {
      ## we'll just use the best to save time
      best_brier_ind <- best_brier_inds[[j]]
      return(permutation_brier(
        x_mat,
        match_list = matches_by_sinks[[j]][[best_brier_ind]]
      ))
    } else {
      return(lapply(matches_by_sinks[[j]], function(match_list) {
        permutation_brier(
          x_mat,
          match_list = match_list
        )
      }))
    }
  })

  ## compute the permutation score for each match
  permutation_brier_scores <- setNames(lapply(seq_len(length(n_sinks)),
                                              function(j) {
                                                if (approximate_by_best) {
                                                  permutation_vec <- permutation_briers[[j]]
                                                  return(unlist(lapply(briers_by_sinks[[j]], function(x) {
                                                    mean(x <= permutation_vec)
                                                  })))
                                                } else {
                                                  return(unlist(lapply(
                                                    seq_len(length(briers_by_sinks[[j]])), function(k) {
                                                      permutation_vec <- permutation_briers[[j]][[k]]
                                                      mean(briers_by_sinks[[j]][[k]] <= permutation_vec)
                                                    })))
                                                }
                                              }), n_sinks)

  ## now that we're doing one-sided brier,
  ## the lowest value will just be the best
  ## so will highest brier

  best_matches <- setNames(lapply(seq_len(length(n_sinks)), function(j) {
    best_brier_ind <- if (approximate_by_best) {
      best_brier_inds[[j]]
    } else {
      which.min(permutation_brier_scores[[j]])
    }

    if (approximate_by_best) {
      ## because they all used the same vector -
      ## so the best must also have the best score here
      stopifnot(
        permutation_brier_scores[[j]][best_brier_ind] ==
          min(permutation_brier_scores[[j]])
      )
    }

    list(
      n_sinks = n_sinks[j],
      raw_brier = briers_by_sinks[[j]][best_brier_ind],
      permutation_brier = permutation_brier_scores[[j]][best_brier_ind],
      match_list = matches_by_sinks[[j]][[best_brier_ind]]
    )
  }), n_sinks)

  list(
    permutation_brier_scores = permutation_brier_scores,
    best_matches = best_matches
  )
}

