% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataGen.R
\name{dataGen}
\alias{dataGen}
\title{Generate Simulation Data}
\usage{
dataGen(N = 10000, n_controls = 10, theta = 0.5, seed_set = NULL)
}
\arguments{
\item{N}{number of desired observations/rows.}

\item{n_controls}{number of control variables/columns.}

\item{theta}{user-defined treatment effect.}

\item{seed_set}{an integer used for set.seed.}
}
\value{
A data frame
}
\description{
This function allows the user to generate a simulation data frame for
evaluating the double machine learning estimator for theta. It requres
the user to specify how many observations and control variables to
include aside from the treatment variable. There is also an
option to set the seed for reproducible results.
}
\examples{
dataGen()
}
