#' swap_pairs
#'
#' Given a logical swap vector, switches treated and control units.
#'
#' @param match_list typical \code{match_list} entity.
#' @param swap logical vector indicating which pairs should swap.
#' @return Another match_list, but now with swapped pairs.
#' @export
swap_pairs <- function(match_list,
                       swap) {
  stopifnot(length(swap) == length(match_list[["treat_index"]]))

  primary_list <- list(
    treat_index = ifelse(!swap,
                         match_list[["treat_index"]],
                         match_list[["control_index"]]
    ),
    control_index = ifelse(swap,
                           match_list[["treat_index"]],
                           match_list[["control_index"]]
    )
  )

  if (all(c("treat_index_within", "control_index_within") %in%
          names(match_list))) {
    primary_list[["treat_index_within"]] <-
      ifelse(!swap,
             match_list[["treat_index_within"]],
             match_list[["control_index_within"]]
      )

    primary_list[["control_index_within"]] <-
      ifelse(swap,
             match_list[["treat_index_within"]],
             match_list[["control_index_within"]]
      )

    primary_list <- primary_list[c(1L, 3L, 2L, 4L)]
  }

  primary_list[["distance"]] <- match_list[["distance"]]

  primary_list
}
