#' Double Machine Learning for Esitmating Treatment Effects
#'
#' This function performs the estimation of the debiased estimator
#' as outlined in Chernozhukov et al (2016). It requres the user to
#' provide a data frame, column indexes for the depenent and treatment
#' variable, how many splits to perform, and whether to summarize
#' the results or return all the output data in a list instead. The user
#' also has a choice of how theta is calculated.
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
dmlRF2 <- function(data    = NULL,    # a data frame
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
  fitMod <- function(j){
    
    # fit the nuisance functions
    g_hat <- randomForest::randomForest(z[setdiff(I_c, j), ],
                                        y[setdiff(I_c, j), ],
                                        maxnodes = ncol(z))
    m_hat <- randomForest::randomForest(z[setdiff(I_c, j), ],
                                        d[setdiff(I_c, j), ],
                                        maxnodes = ncol(z))
    
    # residualize the dep. and treatment vars using nuisance functions
    w <- y[j, ] - predict(g_hat, data.frame(z[j, ]))
    v <- d[j, ] - predict(m_hat, data.frame(z[j, ]))
    
    # compute the treatment coefficient and score function
    if(DML == "DML1"){
      theta <- mean(v * w)/mean(v * d[j])
      score <- (w - (d[j]*theta))*(v)
    } else if(DML == "DML2"){
      theta <- mean(v * w)/mean(v^2)
      score <- (w - (v*theta))*(v)
    } else if(DML == "FWL"){
      theta <- coef(lm(w ~ v + 0))[[1]]
      score <- resid(lm(w ~ v))
    }
    
    # place outcomes into a data frame
    return(data.frame(W = w,
                      V = v,
                      D = d[j],
                      resid = score,
                      Theta = theta,
                      index = j))
  }
  
  # fit the models and get output objects
  lgList = lapply(I, fitMod)
  
  # specify reciever data frame
  out <- data.frame()
  
  # unpack the list into the output data frame
  for(i in 1:length(lgList)){
    out <- rbind(out, lgList[[i]])
  }
  
  # create the summary information as output option
  theta <- mean(out$Theta)
  se    <- sqrt(mean(out$V^2*(out$W - out$V*theta)^2)/mean(out$V^2)^2)/sqrt(N)
  
  if(compile){
    return(c("theta" = theta, "std.err" = se))
  } else {
    return(lgList)
  }
}

