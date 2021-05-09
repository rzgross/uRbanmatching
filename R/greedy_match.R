#' greedy_match
#'
#' Compute greedy matching without replacement. Alternative to "with replacement" and "optimal".
#'
#' @inheritParams bipartite_matches
#' @export
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
