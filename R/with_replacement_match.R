#' with_replacement_match
#'
#' Computes matching with replacement. Alternative to "greedy" and "optimal".
#'
#' @inheritParams bipartite_matches
#' @export
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
