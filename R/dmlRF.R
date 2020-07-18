#' Double Machine Learning for Esitmating Treatment Effects
#'
#' This function performs the estimation of the debiased estimator
#' as outlined in Chernozhukov et al (2016). It requres the user to
#' provide a data frame, column indexes for the depenent and treatment
#' variable, how many splits to perform, and whether to summarize
#' the results or return all the output data in a list instead. The user
#' also has a choice of how theta is calculated.
#'
#' dmlRF differs from dmlRFapply in that the former fits the models to data
#' splits using loops while the latter does so by the apply method.
#'
#' @param data a data frame with depenent, treatment, and control variables.
#' @param dep the column number of the dependent variable.
#' @param treat the column number of the treatment variable.
#' @param compile if TRUE, summarize the results.
#' @param splits number of sample splits used to estimate treatment effect.
#' @param DML 1 - theta as calculated in equation 1.3 of Chernozhukov et al. (2016).
#'            2 - theta as calculated in the footnotes on page 4 of Chernozhukov et al. (2016).
#'            3 - theta using the Frisch-Waugh-Lovell style residual on residual regression.
#' @return A named vector or list
#' @export
#' @examples
#' dmlRF(data, dep = 1, treat = 2, splits = 3, DML = "FWL")

# define dml estimating function -----------------------------------------------
dmlRF <- function(data    = NULL,    # a data frame
                  dep     = NULL,    # col index of dependent variable
                  treat   = NULL,    # col index of treatment variable
                  compile = TRUE,    # summarize results as output?
                  splits  = 2,       # how many sample splits to do
                  DML     = "DML1"){ # which DML estimation to calculate

  # initialize split and receiver objects
  N      <- nrow(data)
  I      <- list()
  I_last <- 1:N
  I_c    <- 1:N
  n_size <- round(N/splits)
  W      <- c()
  V      <- c()
  D      <- c()
  thetas <- c()
  scores <- c()

  # split the dataframe to matrices
  y = as.matrix(data[dep])
  d = as.matrix(data[treat])
  z = as.matrix(data[setdiff(1:ncol(data), c(dep, treat))])

  # perform the splits
  for(i in 1:splits){
    if(i < splits){
      I[[i]] <- sample(I_last, n_size)
      I_last <- setdiff(I_last, I[[i]])
    } else {
      I[[i]] <- I_last
    }
  }

  # iterative estimation of each kth sub-theta
  for(k in 1:splits){
    # fit the nuisance functions
    g_hat <- randomForest::randomForest(z[setdiff(I_c, I[[k]]), ],
                                        y[setdiff(I_c, I[[k]]), ],
                                        maxnodes = ncol(z))
    m_hat <- randomForest::randomForest(z[setdiff(I_c, I[[k]]), ],
                                        d[setdiff(I_c, I[[k]]), ],
                                        maxnodes = ncol(z))

    # residualize the dep. and treatment vars using nuisance functions
    w <- y[I[[k]], ] - predict(g_hat, data.frame(z[I[[k]], ]))
    v <- d[I[[k]], ] - predict(m_hat, data.frame(z[I[[k]], ]))

    # compute the treatment coefficient and score function
    if(DML == "DML1"){
      theta <- mean(v * w)/mean(v * d[I[[k]], ])
      score <- (w - (d[I[[k]], ]*theta))*(v)
    } else if(DML == "DML2"){
      theta <- mean(v * w)/mean(v^2)
      score <- (w - (v*theta))*(v)
    } else if(DML == "FWL"){
      theta <- coef(lm(w ~ v + 0))[[1]]
      score <- resid(lm(w ~ v))
    }

    # running append of output vectors
    W      <- c(W, w)
    V      <- c(V, v)
    D      <- c(D, d[I[[k]], ])
    scores <- c(scores, score)
    thetas <- c(thetas, theta)
  }

  # store all relevant objects in list as output option
  lgList        <- list()
  lgList$splits <- I
  lgList$thetas <- thetas
  lgList$W      <- W
  lgList$V      <- V
  lgList$D      <- D
  lgList$scores <- scores

  # create the summary information as output option
  theta <- mean(thetas)
  se    <- sqrt(mean(V^2*(W - V*theta)^2)/mean(V^2)^2)/sqrt(N)
  t     <- theta/se
  p     <- 2*pt(-abs(t), df = N - 1)

  if(compile){
    return(c("theta" = theta,
             "std.err" = se))
  }else{
    return(lgList)
  }
}
