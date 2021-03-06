% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_caliper_list.R
\name{gen_caliper_list}
\alias{gen_caliper_list}
\title{gen_caliper_list}
\usage{
gen_caliper_list(caliper_vec = NULL, caliper_max = NULL, continuous_mult = 100)
}
\arguments{
\item{caliper_vec}{Default NULL; numeric vector that stops matches if beyond \code{caliper_max}.}

\item{caliper_max}{The maximum allowed difference.}

\item{continuous_mult}{The value to multiply differences above caliper max. Set as \code{Inf} to have infinite penalties.}
}
\value{
Either \code{NULL}, or a list with the same names as the input, after checking values.
}
\description{
Wrapper function to unify caliper input.
}
