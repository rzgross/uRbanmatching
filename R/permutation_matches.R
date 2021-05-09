#' permutation_matches
#'
#' Takes matches and their Brier scores, and computes
#' permutation Brier scores and the best matches.
#'
#' @param matches_by_sinks List by number of sinks, each a list of match results (a match list), for each weight vector.
#' @param briers_by_sinks List by number of sinks, each a vector of Brier results. Basically a number for each match in \code{matches_by_sinks}.
#' @param x_mat Typical input matrix.
#' @param n_sinks Vector of number of sinks.
#' @param approximate_by_best Logical, default \code{TRUE}. Only compute one permutation distribution, using the best result by brier score to do so.
#' @param silent Whether to suppress message output as it's running. Default \code{!interactive()}.
#' @return Returns a list of two lists. The first is vectors of permutation Brier scores (one per match). The second is the best match at each sink value, along with some extra info about that match.
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
