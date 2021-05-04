#' Wrapper to unify tolerance input, along with minor checks
#'
#' Note that we will always block equal equality on the tolerance
#' vec, hence why \code{tolerance_min} being \code{NULL} is equivalent
#' to zero. For more sophisticated behaviour, you may want to
#' directly control the caliper list or weighted distance matrix
#' in the functions that generate matches.
#' @param tolerance_vec Default NULL; numeric "continuous treatment"
#'   vector that we use to form
#'   non-bipartite matches: units i and j can be matched if
#'   \eqn{\mid \texttt{tolerance_vec}[i] - \texttt{tolerance_vec}[j]\mid >
#'   \texttt{tolerance_min}}.
#' @param tolerance_min See above for what this does - blocks
#'   matches that are too close on \code{tolerance_vec}. Something like minimum
#'   relevant difference to be a "treatment" effect. Default \code{NULL} gives
#'   zero, i.e. only blocks equality of \code{tolerance_vec}.
#' @param tolerance_max Optionally we may want to also ensure
#'   our "treatment" values aren't too far apart. E.g. we may think our
#'   assumptions are reasonable valid for small differences
#'   in the tolerance vector, but not for large. Or another way:
#'   we're asking for say marginal effects: how bad is one extra beer a
#'   day?
#' @return Either \code{NULL}, or a list with the same names as the input,
#'   with validated values.
#' @author Colman Humphrey
#'
#' @export
gen_tolerance_list <- function(tolerance_vec = NULL,
                               tolerance_min = NULL,
                               tolerance_max = NULL) {
    if (is.null(tolerance_vec)) {
        if (!is.null(tolerance_min) || !is.null(tolerance_max)) {
            stop(
                "can't set `tolerance_{min, max}` without setting ",
                "`tolerance_vec`"
            )
        }
        return(NULL)
    }

    if (is.null(tolerance_min)) {
        tolerance_min <- 0
    } else {
        if (length(tolerance_min) > 1) {
            stop("`tolerance_min` should be length one (for now) or NULL")
        }

        if (tolerance_min < 0) {
            if (abs(tolerance_min) > sqrt(.Machine$double.eps)) {
                stop("`tolerance_min` should be non-negative")
            }
            ## in this case it's given as negative, but barely
            ## so we'll just use zero
            tolerance_min <- 0
        }
    }

    if (!is.null(tolerance_max)) {
        if (length(tolerance_max) > 1) {
            stop("`tolerance_max` should be length one (for now) or NULL")
        }

        if (tolerance_max < sqrt(.Machine$double.eps)) {
            stop("`tolerance_max` should be strictly positive")
        }

        if (!is.null(tolerance_min) && tolerance_max < tolerance_min) {
            stop(
                "If both given, ",
                "`tolerance_max` must be greater than `tolerance_min`"
            )
        }
    }

    return(list(
        tolerance_vec = tolerance_vec,
        tolerance_min = tolerance_min,
        tolerance_max = tolerance_max
    ))
}


#' Converts tolerance list to caliper list
#' @inheritParams gen_caliper_list
#' @param tolerance_list Result from \code{gen_tolerance_list}
#' @param use_min Logical, should the caliper max be the tolerance min?
#'   Use tolerance max if not. Default \code{TRUE}.
#' @return List from \code{gen_caliper_list}
#' @author Colman Humphrey
#'
#' @keywords internal
tolerance_to_caliper_list <- function(tolerance_list,
                                      use_min = TRUE,
                                      continuous_mult = 1) {
    stopifnot(is_tf(use_min))

    cal_max <- ifelse(use_min,
        tolerance_list[["tolerance_min"]],
        tolerance_list[["tolerance_max"]]
    )
    gen_caliper_list(
        caliper_vec = tolerance_list[["tolerance_vec"]],
        caliper_max = cal_max,
        continuous_mult = continuous_mult
    )
}
