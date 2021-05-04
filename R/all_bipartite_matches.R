#' idea is to complete a search for a given N,
#' or vector N
#'
#' @param x_mat input/design matrix (already rank-adjusted etc)
#' @param cov_x the (potentially rank-adjusted) covariance of \code{x_mat}.
#'   This means it's possible that \code{cov(x_mat)} is not equal to
#'   \code{cov_x}; see \code{covariance_with_ranks} for more details.
#' @param weight_list list of weight vectors. See `generate_random_weights` to
#'   automatically generate a reasonable set of vectors.
#' @param treat_vec Logical (or 1/0) vector, indicating treatment (or control).
#' @param n_sinks how many potential matches to not bother with
#'   NOTE: you can do this as a vector, but not for optimal matching.
#' @param caliper_list Optional, see \code{gen_caliper_list}. Provide
#'   this to force matches that are close on some metric.
#' @param tol_val For optimal matches, you can set a tolerance
#'   to be within optimality of, which can be zero for perfect optimality.
#'   Default 1e-4 is reasonable in many cases.
#'
#' @export
all_bipartite_matches <- function(x_mat,
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
                                  propensity_list = match_propensity_list(NULL),
                                  sqrt_mahal = TRUE,
                                  tol_val = NULL) {
  if (is.null(n_sinks)) {
    n_sinks <- 0L
  }

  if (!is.null(propensity_list)) {
    if (!is.null(caliper_list)) {
      stop(
        "don't use both `caliper_list` and `propensity_list`: ",
        " If you do want both, create the combined caliper separately"
      )
    }

    ## in case of logical
    treat_vec <- treat_vec * 1L

    ## generate propensity score
    prop_list_names <- c(
      "propensity_function",
      "oos_propensity",
      "n_folds"
    )
    prop_score <- propensity_score(
      x_mat = x_mat,
      treat_vec = treat_vec,
      propensity_list = propensity_list[prop_list_names]
    )
    caliper_list <- gen_caliper_list(
      caliper_vec = prop_score,
      caliper_max = sd(prop_score) * propensity_list[["caliper_sd_mult"]],
      continuous_mult = propensity_list[["continuous_mult"]]
    )
  }

  if (!is.null(caliper_list)) {
    caliper_dist_mat <- create_caliper(caliper_list,
                                       treat_vec = treat_vec
    )
  }

  by_weight_list <- lapply(weight_list, function(weight_vec) {
    w_dist_mat <- weighted_mahal(x_mat,
                                 cov_x = cov_x,
                                 weight_vec = weight_vec,
                                 treat_vec = treat_vec,
                                 sqrt_mahal = sqrt_mahal
    )

    if (!is.null(caliper_list)) {
      w_dist_mat <- w_dist_mat + caliper_dist_mat
    }

    bipartite_matches(
      dist_mat = w_dist_mat,
      treat_vec = treat_vec,
      match_method = match_method,
      n_sinks = n_sinks,
      tol_val = tol_val,
      weight_vec = weight_vec
    )
  })

  ## more natural to group by sink value
  setNames(lapply(n_sinks, function(x) {
    lapply(by_weight_list, function(y) {
      y[[as.character(x)]]
    })
  }), n_sinks)
}



#' Generating bipartite matched pairs
#'
#' Generates matched pairs either:
#' \describe{
#'   \item{With Replacement}{Finds smallest control for each treatment}
#'   \item{Without Replacement, Greedy}{Greedily generates pairs. Note that
#'   the order for choosing the greedy pairs is random, which is not the only
#'   possible solution.}
#'   \item{Without Replacement, Optimally}{Minimum total distance}
#' }
#' If you're happy to use control units potentially multiple times,
#' then the first way is fast and optimal.
#'
#' If not, you have to trade off speed vs optimality. Greedy runs
#' over all units in a random order, so if you want to run greedy a bunch of
#' times and take the best, it would still be (likely) much faster than
#' running optimal matching.
#'
#' @param dist_mat Matrix of pairwise distances.
#' @param treat_vec Vector representing all subjects; 0 for control,
#'   1 for treated.
#' @param match_method This enum corresponds to the three matching methods
#'   discussed above:
#'   \describe{
#'     \item{"with_replacement"}{Finds smallest control for each treatment}
#'     \item{"greedy"}{Greedily generates pairs. Note that
#'     the order for choosing the greedy pairs is random, which is not the only
#'     possible solution.}
#'     \item{"optimal"}{Minimum total distance}
#'   }
#' @param n_sinks how many sinks to use; can be vector. Note that for greedy and
#'   simple with-replacement matching, it's often better to sort this
#'   elsewhere. Optimal matching can only take one value.
#'   Default NULL to match all and ignore sinks.
#' @param tol_val tolerance for solving optimal matches - how far is
#'   is acceptable to be from the true optimal value? Speed with large value,
#'   accuracy with small. Only relevant for \code{!with_replacement && !greedy}.
#'   Default 1e-4 is reasonable in many cases.
#' @param weight_vec Default \code{NULL}: optionally supply the weight vector
#'   used to generate \code{dist_mat} and it'll be returned in the
#'   \code{match_list} generated from this function
#' @return basic return value is a list with five elements and an optional sixth:
#'   \describe{
#'     \item{\code{treat_index}}{index of treated units, from all units}
#'     \item{\code{treat_index_within}}{index of treated units,
#'           from the set of treated}
#'     \item{\code{control_index}}{index of control units, from all units}
#'     \item{\code{control_index_within}}{index of control units,
#'           from the set of control}
#'     \item{\code{distance}}{distances between the pairs}
#'     \item{\code{weight_vec}}{weight vector used to generate
#'           \code{dist_mat} if supplied}
#'   }
#'   You'll get a list of such objects, each
#'   with an extra element: the number of sinks used. If you
#'   have \code{n_sinks} as \code{NULL}, then it'll default to
#'   a single sink value: zero.
#'
#' @export
bipartite_matches <- function(dist_mat,
                              treat_vec,
                              match_method = c(
                                "with_replacement",
                                "optimal",
                                "greedy"
                              ),
                              n_sinks = NULL,
                              tol_val = NULL,
                              weight_vec = NULL) {
  stopifnot(is.matrix(dist_mat))
  stopifnot(min(dist_mat) >= 0)

  match_method <- match.arg(match_method)

  ## in case of logical, or double
  treat_vec <- as.integer(treat_vec * 1L)

  stopifnot(all(unique(treat_vec) %in% c(0L, 1L)))
  stopifnot(length(treat_vec) == (nrow(dist_mat) + ncol(dist_mat)))

  if (!is.null(n_sinks)) {
    stopifnot(is.numeric(n_sinks) &&
                min(n_sinks) >= 0L &&
                !any(is.na(n_sinks)) &&
                length(unique(n_sinks)) == length(n_sinks))
  }

  if (match_method != "optimal") {
    if (!is.null(tol_val)) {
      stop("tol_val should only be set for optimal matching",
           call. = FALSE
      )
    }
  } else {
    tol_val <- ifelse(!is.null(tol_val),
                      tol_val, 1e-4
    )
    stopifnot(is.numeric(tol_val) && length(tol_val) == 1L &&
                !is.na(tol_val))
  }

  ## ------------------------------------

  if (match_method == "with_replacement") {
    return(simple_sink_wrap(
      with_replacement_match(
        dist_mat,
        treat_vec
      ),
      n_sinks,
      weight_vec
    ))
  }

  if (match_method == "greedy") {
    return(simple_sink_wrap(
      greedy_match(
        dist_mat,
        treat_vec
      ),
      n_sinks,
      weight_vec
    ))
  }

  optimal_sink_wrap(
    dist_mat,
    treat_vec,
    n_sinks,
    tol_val,
    weight_vec
  )
}


#' @inheritParams bipartite_matches
#' @keywords internal
with_replacement_match <- function(dist_mat,
                                   treat_vec) {
  control_index <- which(treat_vec == 0L)

  match_list <- list(
    treat_index = which(treat_vec == 1L),
    treat_index_within = seq_len(sum(treat_vec))
  )

  matched_vec <- apply(dist_mat, 1, function(x) {
    which(rank(x, ties.method = "random") == 1L)
  })

  match_list[["control_index"]] <- control_index[matched_vec]
  match_list[["control_index_within"]] <- matched_vec
  match_list[["distance"]] <- dist_mat[cbind(
    seq_len(nrow(dist_mat)),
    matched_vec
  )]

  match_list
}


#' @inheritParams bipartite_matches
#' @keywords internal
greedy_match <- function(dist_mat,
                         treat_vec) {
  control_index <- which(treat_vec == 0L)

  match_list <- list(
    treat_index = which(treat_vec == 1L),
    treat_index_within = 1:sum(treat_vec)
  )

  min_vals <- apply(dist_mat, 1, min)

  result_mat <- matrix(NA, nrow = nrow(dist_mat), ncol = 3)

  while (min(min_vals) < Inf) {
    random_value <- sample(1:length(min_vals),
                           size = 1,
                           prob = 1 / (min_vals + 1)
    )
    match_ind <- which(rank(dist_mat[random_value, ],
                            ties.method = "random"
    ) == 1L)
    result_mat[random_value, 3] <- dist_mat[random_value, match_ind]
    result_mat[random_value, 1] <- match_ind
    result_mat[random_value, 2] <- control_index[match_ind]

    ## blocked:
    dist_mat[, match_ind] <- Inf
    dist_mat[random_value, ] <- Inf
    min_vals <- apply(dist_mat, 1, min)
  }
  result_mat <- result_mat[!is.na(result_mat[, 1]), ]

  match_list[["control_index_within"]] <- result_mat[, 1, drop = TRUE]
  match_list[["control_index"]] <- result_mat[, 2, drop = TRUE]
  match_list[["distance"]] <- result_mat[, 3, drop = TRUE]

  match_list
}


#' Computes optimal matches, bipartite
#'
#' @inheritParams bipartite_matches
#' @param n_sinks single value
#' @keywords internal
optimal_match <- function(dist_mat,
                          treat_vec,
                          n_sinks = 0,
                          tol_val = 1e-4) {
  ## ------------------------------------

  dist_mat_sink <- cbind(
    dist_mat,
    ## sinks are zeros
    matrix(0, nrow = nrow(dist_mat), ncol = n_sinks)
  )
  main_treat_index <- which(treat_vec == 1L)
  rownames(dist_mat_sink) <- main_treat_index

  all_treated <- treat_vec
  if (n_sinks > 0) {
    ## the sinks are not treated
    all_treated <- c(treat_vec, rep(0L, n_sinks))
  }
  all_control_index <- which(all_treated == 0L)
  colnames(dist_mat_sink) <- all_control_index

  treat_data <- data.frame(treated = all_treated)
  rownames(treat_data) <- seq_len(nrow(treat_data))

  ## Matching using optmatch
  ## if need be:
  ## options("optmatch_max_problem_size" = Inf))
  match_vec <- optmatch::pairmatch(dist_mat_sink,
                                   tol = tol_val,
                                   data = treat_data
  )
  list_results <- aggregate(names(match_vec),
                            by = list(match_vec),
                            FUN = function(x) {
                              as.numeric(x)
                            },
                            simplify = FALSE
  )[["x"]]
  first_vec <- unlist(lapply(list_results, function(x) x[1]))
  second_vec <- unlist(lapply(list_results, function(x) x[2]))
  control_ind <- all_treated[first_vec] == 0
  keep_ind <- pmax(first_vec, second_vec) <= length(treat_vec)

  treat_subj_ind <- ifelse(control_ind, second_vec, first_vec)[keep_ind]
  control_subj_ind <- ifelse(control_ind, first_vec, second_vec)[keep_ind]

  ord_vec <- order(treat_subj_ind)

  treat_subj_ind <- treat_subj_ind[ord_vec]
  control_subj_ind <- control_subj_ind[ord_vec]

  match_list <- list(
    treat_index = treat_subj_ind,
    treat_index_within = match(treat_subj_ind, main_treat_index),
    control_index = control_subj_ind,
    control_index_within = match(control_subj_ind, all_control_index)
  )
  match_list[["distance"]] <- dist_mat[cbind(
    match_list[["treat_index_within"]],
    match_list[["control_index_within"]]
  )]

  match_list
}


#' Simple wrapper to unify caliper input
#'
#' @param caliper_vec Default NULL; numeric vector that "blocks"
#'   matches if they're too (further than \code{caliper_max}) on this value
#' @param caliper_max The maximum allowed difference (exactly this difference
#'   is allowed).
#' @param continuous_mult The value to multiply differences above caliper max.
#'   Set as \code{Inf} to have infinite penalties, i.e. block matches above
#'   the max.
#' @return Either \code{NULL}, or a list with the same names as the input, after checking
#'   values.
#'
#' @export
gen_caliper_list <- function(caliper_vec = NULL,
                             caliper_max = NULL,
                             continuous_mult = 100) {
  if (is.null(caliper_vec)) {
    if (!is.null(caliper_max)) {
      stop("can't give `caliper_max` without `caliper_vec`")
    }
    return(NULL)
  }

  if (is.null(caliper_max)) {
    stop("supply `caliper_max` if using calipers")
  }

  if (length(caliper_max) > 1L) {
    stop("`caliper_max` should be length one")
  }

  if (length(continuous_mult) > 1L) {
    stop("`continuous_mult` should be length one")
  }

  return(list(
    caliper_vec = caliper_vec,
    caliper_max = caliper_max,
    continuous_mult = continuous_mult
  ))
}


#' Creates caliper penalties
#'
#' Given a vector, \code{caliper_vec}, this function sets penalties
#' for any pairwise differences above \code{caliper_max}; either
#' an \code{Inf} penalty (the default), or a continuous penalty. Typical usage will
#' then add the resulting matrix from this function onto a distance
#' matrix, say from pairwise Mahalanobis.
#'
#' @param caliper_list Result of \code{gen_caliper_list}
#' @param treat_vec Optional; if you only want pairs between treat and control.
#' @return A matrix, either square (\code{|caliper_vec| x |caliper_vec|})
#'   or else \code{sum(treat_vec == 1) x sum(treat_vec == 0)}.
#'
#' @export
create_caliper <- function(caliper_list,
                           treat_vec = NULL) {
  if (is.null(treat_vec)) {
    caliper_treat <- caliper_list[["caliper_vec"]]
    caliper_control <- caliper_list[["caliper_vec"]]
  } else {
    if (!(length(treat_vec) == length(caliper_list[["caliper_vec"]]))) {
      stop("treat_vec not the same length as `caliper_vec`")
    }
    ## in case it's logical...:
    treat_vec <- treat_vec * 1L

    caliper_treat <- caliper_list[["caliper_vec"]][treat_vec == 1]
    caliper_control <- caliper_list[["caliper_vec"]][treat_vec == 0]
  }

  abs_caliper_diff <- abs(outer(caliper_treat,
                                caliper_control,
                                FUN = "-"
  )) - caliper_list[["caliper_max"]]
  abs_caliper_diff[abs_caliper_diff < 0] <- 0

  continuous_mult <- caliper_list[["continuous_mult"]]

  if (!is.null(continuous_mult) && continuous_mult < Inf) {
    penalty_mat <- continuous_mult * abs_caliper_diff
  } else {
    penalty_mat <- abs_caliper_diff
    penalty_mat[penalty_mat > 0] <- Inf
  }

  penalty_mat
}



#' Generates the propensity parameters used for using propensity-based calipers
#'
#' We use this for input for \code{all_propensity_caliper_matches} (and
#' likely a nonbipartite version soon).
#' @inheritParams gen_propensity_list
#' @param caliper_sd_mult We'll set the maximum gap between units
#'   as \code{sd(propensity_score) * k}, where this parameter is the
#'   value k. Default 0.6.
#' @param continuous_mult See e.g. \code{gen_caliper_list}: instead of
#'   blocking matches that are "too far apart" on the caliper, we'll
#'   add a penalty for going above.
#' @return list with names equal to all input params
#' @export
match_propensity_list <- function(propensity_function = propensity_score_xgb(),
                                  oos_propensity = FALSE,
                                  n_folds = NULL,
                                  caliper_sd_mult = 0.6,
                                  continuous_mult = 100) {
  if (is.null(propensity_function)) {
    return(NULL)
  }

  plain_prop_list <- gen_propensity_list(
    propensity_function = propensity_function,
    oos_propensity = oos_propensity,
    n_folds = n_folds
  )

  list(
    propensity_function = plain_prop_list[["propensity_function"]],
    oos_propensity = plain_prop_list[["oos_propensity"]],
    n_folds = plain_prop_list[["n_folds"]],
    caliper_sd_mult = caliper_sd_mult,
    continuous_mult = continuous_mult
  )
}


#' Calculates propensity scores for a given matrix and treatment vector
#'
#' This function takes in an input matrix and a treatment vector,
#' along with a function that makes predictions (default xgboost method given)
#' and returns a predicted probability of treatment for each unit, either
#' using in-sample or out-of-sample fits.
#' @param x_mat Standard input matrix (already rank adjusted).
#' @param treat_vec Usual 0/1 treatment vector.
#' @param propensity_list See \code{gen_propensity_list}
#' @return Returns a vector equal in length to treat_vec of propensity.
#'   score.
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


#' Computes weighted Mahalanobis distance, using Choleski decomp.
#'
#' Note: def. of weighted:
#' \eqn{d(x_i, x_j) = (x_i - x_j)' W \Sigma^{-1} W (x_i - x_j)}
#' where \eqn{W} = \code{diag(weight_vec)}
#' R's cholesky gives \eqn{U} s.t. \eqn{U' U = S} (i.e. \eqn{U = chol(S)})
#' so in general we want:
#' \deqn{
#'   (x_i - x_j)' W (U' U)^{-1} W (x_i - x_j) \\
#' = (x_i - x_j)' W U^{-1} (U')^{-1} W (x_i - x_j) \\
#' = x_i' W U^{-1} (U')^{-1} W x_i +
#'   x_j' W U^{-1} (U')^{-1} W x_j \\
#'   - 2 x_i' W U^{-1} (U')^{-1} W x_j
#' }
#' Solving the above is manageable if we have \eqn{y_i = (U')^{-1} W x_i}
#' or moving terms around, \eqn{W^{-1} U' y_i = x_i},
#' which is simple by \code{forwardsolve}. Then we have:
#' \deqn{dist(x_i, x_j) = ||y_i||^2 + ||y_j||^2 - 2 y_i' y_j}
#' Letting \eqn{Y = (y_1', y_2', .... ,y_n')'}
#' (which we can get in one line,
#' \eqn{Y' = } \code{forwardsolve} \eqn{(W^{-1} U', x_mat')})
#' the first two parts are just \code{rowSums}\eqn{(Y^2)}
#' (note that the code uses \eqn{Y'}, thus \code{colSums})
#' (add \code{outer(.,.)} to finish)
#' and the last is \eqn{Y Y'}
#'
#' @param x_mat numeric matrix (adjust non-numeric columns prior),
#'   already rank-adjusted if desired
#' @param cov_x covariance of x, calculated potentially with ranks
#' @param weight_vec vector of weights corresponding to columns of x_mat,
#'   giving weights relative to "raw" (ranked) Mahalanobis.
#'   Note that the resulting matrix does depend on the scale of the
#'   weight vector, but the matching won't: scaling the Mahalanobis
#'   matrix has no effect on distance minimising pairs etc
#' @param treat_vec optionally specify which units are treated.
#'   If \code{NULL} (default), will just return
#'   \code{nrow(x_mat) x nrow(x_mat)} distance
#'   matrix of all pairs. Can be logicals or \eqn{{0, 1}}
#' @param sqrt_mahal logical, default TRUE; do you want regular Mahalanobis:
#'   \eqn{d(x_i, x_j) = (x_i - x_j)' \Sigma^{-1} (x_i - x_j)}
#'   or the square root? (in weighted, \eqn{\Sigma^{-1}} becomes
#'   \eqn{W \Sigma^{-1} W} in the above)
#' @param partial_index In some cases, you want a subset of the full
#'   N x N matrix, but you can't partition it like you want
#'   with treat_vec. e.g. you want
#'   \code{full_dist[c(1, 2, 3), 1:10]}
#'   then use
#'   \code{partial_index = list(c(1,2,3), 1:10)}
#' @return returns a matrix of pairwise distances; the relevant indexing
#'   depends on \code{treat_vec} and \code{partial_index}
#' @export
weighted_mahal <- function(x_mat,
                           cov_x,
                           weight_vec = NULL,
                           treat_vec = NULL,
                           sqrt_mahal = TRUE,
                           partial_index = NULL) {
  if (!is.null(partial_index) && !is.null(treat_vec)) {
    stop("Supply at most one of treat_vec and partial_index", call. = FALSE)
  }
  if (!is.matrix(x_mat)) {
    stop("x_mat must be a matrix", call. = FALSE)
  }

  chol_cov <- chol(cov_x)

  if (is.null(weight_vec)) {
    weight_vec <- rep(1 / ncol(x_mat), times = ncol(x_mat))
  } else {
    if (length(weight_vec) != ncol(x_mat)) {
      stop("`weight_vec` should have length equal to `ncol(x_mat)`")
    }
  }

  weighted_ymat <- forwardsolve(diag(1 / weight_vec) %*%
                                  t(chol_cov), t(x_mat))

  ymat_sumsq <- colSums(weighted_ymat^2)

  if (!is.null(treat_vec) || !is.null(partial_index)) {
    if (!is.null(treat_vec)) {
      ## in case it's logical...:
      treat_vec <- treat_vec * 1

      row_index <- treat_vec == 1
      col_index <- treat_vec == 0
    } else {
      row_index <- partial_index[[1]]
      col_index <- partial_index[[2]]
    }

    crossprod_weighted <- t(weighted_ymat[, row_index]) %*%
      weighted_ymat[, col_index]

    mahal_mat <- outer(ymat_sumsq[row_index],
                       ymat_sumsq[col_index],
                       FUN = "+"
    ) -
      2 * crossprod_weighted
  } else {
    crossprod_weighted <- t(weighted_ymat) %*% weighted_ymat

    mahal_mat <- outer(ymat_sumsq, ymat_sumsq, FUN = "+") -
      2 * crossprod_weighted
  }

  ## correcting for e.g. system tol (well, working with doubles)
  mahal_mat <- pmax(mahal_mat, 0)

  if (sqrt_mahal) {
    mahal_mat <- sqrt(mahal_mat)
  }

  mahal_mat
}


#' Wraps a match_list from a simple match result and slices off
#' various numbers of sinks
#'
#' Note that for \code{with_replacement} matching, this is the
#' right thing to do, but for greedy matching with a random ordering
#' it isn't obviously correct, but does give the nice property
#' of more sinks = lower distances.
#' @param simple_match_list match result from one of \code{with_replacement_match},
#'   \code{greedy_match}, \code{with_replacement_nbp_match} or
#'   \code{greedy_nbp_match}.
#' @param n_sinks default NULL, vector of sink values to use.
#' @param weight_vec Default \code{NULL}: optionally supply the weight vector
#'   used to generate \code{dist_mat} and it'll be returned in the
#'   \code{match_list} generated from this function
#' @return list of lists; see parent function
#'
#' @keywords internal
simple_sink_wrap <- function(simple_match_list,
                             n_sinks = NULL,
                             weight_vec = NULL) {
  if (is.null(n_sinks)) {
    n_sinks <- 0L
  }

  ## ------------------------------------

  dist_ranks <- rank(simple_match_list[["distance"]],
                     ties.method = "random"
  )
  setNames(lapply(n_sinks, function(sink_val) {
    keep_ind <- dist_ranks <= (length(dist_ranks) - sink_val)
    match_list <- lapply(simple_match_list, function(x) {
      x[keep_ind]
    })
    match_list[["num_sinks"]] <- sink_val
    if (!is.null(weight_vec)) {
      match_list[["weight_vec"]] <- weight_vec
    }
    match_list
  }), n_sinks)
}


#' Given a vector of sink values, generates an optimal match
#' for each.
#'
#' Will be slow; you can't just generate one match and subset from it.
#' @inheritParams bipartite_matches
#' @param n_sinks default NULL, vector of sink values to use.
#' @param weight_vec Default \code{NULL}: optionally supply the weight vector
#'   used to generate \code{dist_mat} and it'll be returned in the
#'   \code{match_list} generated from this function
#' @return list of lists; see parent function
#'
#' @keywords internal
optimal_sink_wrap <- function(dist_mat,
                              treat_vec,
                              n_sinks,
                              tol_val,
                              weight_vec = NULL) {
  if (is.null(n_sinks)) {
    n_sinks <- 0L
  }

  ## ------------------------------------

  setNames(lapply(n_sinks, function(sink_val) {
    match_list <- optimal_match(
      dist_mat,
      treat_vec,
      sink_val,
      tol_val
    )
    match_list[["num_sinks"]] <- sink_val
    if (!is.null(weight_vec)) {
      match_list[["weight_vec"]] <- weight_vec
    }
    match_list
  }), n_sinks)
}
