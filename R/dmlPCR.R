#' Double Machine Learning for Estimating Treatment Effects
#'
#' One of several functions in this package that performs the estimation of the
#' debiased estimator as outlined in Chernozhukov et al (2016). It requires the
#' user to provide a data frame, column indexes for the dependent and treatment
#' variable, how many splits to perform, and whether to summarize the results or
#' return all the output data in a list instead. The user also has a choice of
#' how theta is calculated.
#'
#' dmlPCR specifically does this using principal component regression from the
#' pls package as the foundational ML model for fitting the nuisance function in
#' cross-splitting.
#'
#' @param data a data frame with dependent, treatment, and control variables.
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
#' dmlPCR(data, dep = 1, treat = 2, splits = 3, DML = "FWL")

# define dml estimating function -----------------------------------------------
dmlPCR <- function(data    = NULL,
                   dep     = NULL,
                   treat   = NULL,
                   compile = TRUE,
                   splits  = 2,
                   DML     = "DML1"){
  
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
  
  # perform the splits
  for(i in 1:splits){
    if(i < splits){
      I[[i]] <- sample(I_last, n_size)
      I_last <- setdiff(I_last, I[[i]])
    } else {
      I[[i]] <- I_last
    }
  }
  
  # adjust the column names for easier downstream processing
  colnames(data)[dep]   = "y"
  colnames(data)[treat] = "d"
  
  # iterative estimation of each kth sub-theta
  for(k in 1:splits){
    # specify the formulas for the nuicanse functions
    gform = paste("y ~",
                  paste(colnames(data)[-c(dep, treat)], collapse = " + "),
                  collapse = " ")
    mform = paste("d ~",
                  paste(colnames(data)[-c(dep, treat)], collapse = " + "),
                  collapse = " ")
    
    # coss-validate and fit the principal componenet regressions
    ghat = pls::pcr(as.formula(gform),
                    data = data[setdiff(I_c, I[[k]]), ],
                    scale = TRUE,
                    validation = "CV")
    mhat = pls::pcr(as.formula(mform),
                    data = data[setdiff(I_c, I[[k]]), ],
                    scale = TRUE,
                    validation = "CV")
    
    # get the ideal number of components
    g_comp = pls::selectNcomp(ghat, method = "onesigma")
    m_comp = pls::selectNcomp(mhat, method = "onesigma")
    
    # get the predicted values
    w <- data[I[[k]], "y"] - predict(ghat,
                                     data.frame(data[I[[k]], -c(dep, treat)]),
                                     ncomp = g_comp)
    v <- data[I[[k]], "d"] - predict(mhat,
                                     data.frame(data[I[[k]], -c(dep, treat)]),
                                     ncomp = m_comp)
    
    # compute the treatment coefficient and score function
    if(DML == "DML1"){
      theta <- mean(v * w)/mean(v * data[I[[k]], "d"])
      score <- (w - (data[I[[k]], "d"]*theta))*(v)
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
    D      <- c(D, data[I[[k]], "d"])
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



