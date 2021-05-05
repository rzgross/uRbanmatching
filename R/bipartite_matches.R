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
