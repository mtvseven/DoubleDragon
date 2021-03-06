% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmlRF.R
\name{dmlRF}
\alias{dmlRF}
\title{Double Machine Learning for Estimating Treatment Effects}
\usage{
dmlRF(
  data = NULL,
  dep = NULL,
  treat = NULL,
  compile = TRUE,
  splits = 2,
  DML = "DML1"
)
}
\arguments{
\item{data}{a data frame with dependent, treatment, and control variables.}

\item{dep}{the column number of the dependent variable.}

\item{treat}{the column number of the treatment variable.}

\item{compile}{if TRUE, summarize the results.}

\item{splits}{number of sample splits used to estimate treatment effect.}

\item{DML}{1 - theta as calculated in equation 1.3 of Chernozhukov et al. (2016).
2 - theta as calculated in the footnotes on page 4 of Chernozhukov et al. (2016).
3 - theta using the Frisch-Waugh-Lovell style residual on residual regression.}
}
\value{
A named vector or list
}
\description{
One of several functions in this package that performs the estimation of the
debiased estimator as outlined in Chernozhukov et al (2016). It requires the
user to provide a data frame, column indexes for the dependent and treatment
variable, how many splits to perform, and whether to summarize the results or
return all the output data in a list instead. The user also has a choice of
how theta is calculated.
}
\details{
dmlRF specifically does this using the random forest from the randomForest
package as the foundational ML model for fitting the nuisance function in
cross-splitting.
}
\examples{
dmlRF(data, dep = 1, treat = 2, splits = 3, DML = "FWL")
}
