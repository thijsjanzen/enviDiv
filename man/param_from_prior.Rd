% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/param_from_prior.R
\name{param_from_prior}
\alias{param_from_prior}
\title{draw parameter combinations from the prior}
\usage{
param_from_prior(empty_input)
}
\arguments{
\item{empty_input}{empty vector, used to facilitate usage of apply}
}
\value{
vector with 7 entries: extinction, sympatric speciation at high
water, sympatric speciation at low water, allopatric speciation, amount of
perturbation, the chosen water model and the abc-weight
}
\description{
draw parameter combinations from the prior
}
