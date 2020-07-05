#' Generate Simulation Data
#'
#' This function allows the user to generate a simulation data frame for
#' evaluating the double machine learning estimator for theta. It requres
#' the user to specify how many observations and control variables to
#' include aside from the treatment variable. There is also an
#' option to set the seed for reproducible results.
#'
#' @param N number of desired observations/rows.
#' @param n_controls number of control variables/columns.
#' @param theta user-defined treatment effect.
#' @param seed_set an integer used for set.seed.
#' @return A data frame
#' @export
#' @examples
#' dataGen()

dataGen <- function(N = 10000,        # how many observations to create
                    n_controls = 10,  # how many controls to create
                    theta = 0.5,      # specify the treatment effect
                    seed_set = NULL){ # set the seed for replicability

  # if applicable, set the seed
  set.seed(seed_set)

  # Generate covariance matrix of z
  k = n_controls
  b = 1/(1:k)
  sigma = clusterGeneration::genPositiveDefMat(k,"unifcorrmat")$Sigma
  sigma = cov2cor(sigma)

  # generate all the data components
  z = mvtnorm::rmvnorm(N, sigma = sigma)
  g = as.vector(cos( z %*% b )^2)
  m = as.vector(sin(z %*% b) + cos(z %*% b))
  d = m + rnorm(N)
  y = theta * d + g + rnorm(N)

  # return the data
  return(data.frame(y, d, z))

}
