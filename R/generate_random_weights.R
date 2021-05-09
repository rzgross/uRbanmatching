#' generate_random_weights
#'
#' Generating random weights to use in matching procedures,
#' each summing to one.
#'
#' @param prior_weights Must be equal to the length of your columns,
#'   i.e. the length of the weight vectors this function will produce.
#' @param number_vectors How many weight vectors to generate.
#' @param minimum_weights If you want to set minimums weights.
#' @param hierarchical_list List per group / category of variable:
#'   \describe{
#'     \item{\code{"index"}}{Vector of indices that this group corresponds with}
#'     \item{\code{"weight"}}{Weight for this group}
#'     \item{\code{"variance"}}{(Optional) how much variance this group will havve}
#'   }
#' @return list of weight vectors
#'
#' @export
generate_random_weights <- function(prior_weights,
                                    number_vectors,
                                    minimum_weights = NULL,
                                    hierarchical_list = NULL) {
  if (is.null(minimum_weights)) {
    minimum_weights <- rep(0, length(prior_weights))
  } else {
    if (length(minimum_weights) == 1) {
      minimum_weights <- rep(minimum_weights, length(prior_weights))
    }

    if (length(minimum_weights) != length(prior_weights)) {
      stop(
        "`minimum_weights` must be length 1 ",
        "or same length as `prior_weights`"
      )
    }
  }

  if (sum(minimum_weights) > 1) {
    stop(
      "`sum(minimum_weights)` must be less than 1, ",
      "since we're returning scaled weights"
    )
  }

  stopifnot(length(number_vectors) == 1L &&
              is.numeric(number_vectors) &&
              number_vectors >= 1L)

  ## ------------------------------------

  if (!is.null(hierarchical_list)) {
    all_index <- unlist(lapply(hierarchical_list, function(x) {
      x[["index"]]
    }))

    if (!isTRUE(all.equal(
      sort(all_index),
      1L:length(prior_weights)
    ))) {
      stop(
        "the hierarchical list must exactly cover ",
        "`1L:length(prior_weights)`"
      )
    }

    hierarchical_list <- lapply(hierarchical_list, function(x) {
      if (is.null(x[["variance"]])) {
        ## scale invariance
        x[["gamma_alpha"]] <- 1
        x[["gamma_beta"]] <- 1 / x[["weight"]]
      } else {
        x[["gamma_beta"]] <- x[["weight"]] / x[["variance"]]
        x[["gamma_alpha"]] <- x[["weight"]] * x[["gamma_beta"]]
      }
      return(x)
    })

    return(hierarchical_random_weights(
      prior_weights = prior_weights,
      number_vectors = number_vectors,
      minimum_weights = minimum_weights,
      hierarchical_list
    ))
  }

  rescale_to_min <- (function(sum_weights, minimum_weights) {
    function(weight_vec) {
      (weight_vec / sum(weight_vec)) *
        (1 - sum_weights) + minimum_weights
    }
  })(sum(minimum_weights), minimum_weights)

  lapply(1L:number_vectors, function(j) {
    temp_weight <- runif(length(prior_weights), 0, prior_weights)
    rescale_to_min(temp_weight)
  })
}
