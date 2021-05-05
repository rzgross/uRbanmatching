#' Constructs an \code{index_list} from a \code{match_list} and a logical
#' index for training data
#'
#' @param match_list see \code{bipartite_matches} etc
#' @param train_index logical index, same length as
#'   \code{match_list[["treat_index"]]} (and
#'   \code{match_list[["control_index"]]})
#' @return returns an \code{index_list} object, see
#'   \code{predict_prepare}
#'
#' @export
index_list_from_match <- function(match_list,
                                  train_index) {
  list(
    treat_train = match_list[["treat_index"]][train_index],
    control_train = match_list[["control_index"]][train_index],
    treat_test = match_list[["treat_index"]][!train_index],
    control_test = match_list[["control_index"]][!train_index]
  )
}
