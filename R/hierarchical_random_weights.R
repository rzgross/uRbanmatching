#' hierarchical_random_weights
#'
#' Generates weights in a hierarchical setting
#'
#' @inheritParams generate_random_weights
#' @export
hierarchical_random_weights <- function(prior_weights,
                                        number_vectors,
                                        minimum_weights,
                                        hierarchical_list) {
  index_list <- lapply(hierarchical_list, function(cat) {
    cat[["index"]]
  })

  order_vec <- order(unlist(index_list))

  rescale_to_min <- (function(sum_weights, minimum_weights) {
    function(weight_vec) {
      (weight_vec / sum(weight_vec)) *
        (1 - sum_weights) + minimum_weights
    }
  })(sum(minimum_weights), minimum_weights)

  lapply(1L:number_vectors, function(j) {
    weight_list <- lapply(hierarchical_list, function(cat) {
      gamma_var <- rgamma(1,
                          shape = cat[["gamma_alpha"]],
                          rate = cat[["gamma_beta"]]
      )
      pre_weights <- runif(
        length(cat[["index"]]),
        0,
        prior_weights[cat[["index"]]]
      )
      (pre_weights / sum(pre_weights)) * gamma_var
    })
    weight_vec <- unlist(weight_list)[order_vec]

    rescale_to_min(weight_vec)
  })
}
