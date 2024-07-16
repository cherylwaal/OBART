#' obart
#'
#' @param x.train Explanatory variables for training data.
#' @param y.train Binary dependent variable for training (in sample) data.
#' @param x.test Explanatory variables for test data.
#' @param xinfo A list (matrix) where the items (rows) are the covariates and 
#' the contents of the items (columns) are the cutpoints.
#' @param usequants If usequants = FALSE, then the cutpoints in xinfo are generated uniformly; 
#' otherwise, if TRUE, uniform quantiles are used for the cutpoints.
#' @param cont Whether or not to assume all variables are continuous.
#' @param rm.const Whether or not to remove constant variables.
#' @param k The number of prior standard deviations \eqn{f(x)} is away from \eqn{+/-3}. 
#' The bigger k is, the more conservative the fitting will be.
#' @param power Power parameter for tree prior.
#' @param base Base parameter for tree prior.
#' @param rho Sparse parameter: typically \eqn{\rho=p} where \eqn{p} is the number of 
#' covariates under consideration.
#' @param ntree The number of trees in the sum.
#' @param numcut The number of possible values of c. 
#' @param ndpost The number of posterior draws returned.
#' @param nskip Number of MCMC iterations to be treated as burn in.
#' @param oldgamma Initial value of cutoff values.
#' @param oldz Initial value of latent continuous variables.
#' @param printevery As the MCMC runs, a message is printed every printevery draws.
#' 
#' @export

obart=function(
    x.train, y.train, x.test=matrix(0.0,0,0),
    xinfo=matrix(0.0,0,0), usequants=FALSE,
    cont=FALSE, rm.const=TRUE, 
    k=2.0, power=2.0, base=.95,rho=NULL,
    ntree=200L, numcut=100L,
    ndpost=1000L, nskip=100L,oldgamma,oldz,
    printevery=100,offset
)
{
  #--------------------------------------------------
  #data
  n = length(y.train)
  
  ## if(binaryOffset!=0)
  ##     stop('binaryOffset not supported by lbart')
  
    temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                           cont=cont, xinfo=xinfo, rm.const=rm.const)
    x.train = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    ## if(length(x.test)>0)
    ##         x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
    if(length(x.test)>0) {
      x.test = bartModelMatrix(x.test)
      x.test = t(x.test[ , temp$rm.const])
    }
    rm.const <- temp$rm.const
    rm(temp)
  
  if(n!=ncol(x.train))
    stop('The length of y.train and the number of rows in x.train must be identical')
  
  p = nrow(x.train)
  np = ncol(x.test)
  if(length(rho)==0) rho <- p
  if(length(rm.const)==0) rm.const <- 1:p
  
  tau=3/(k*sqrt(ntree))

  #--------------------------------------------------
  #call
  res = .Call("cobart", 
               n,#number of observations in training data
               p,  #dimension of x
               np, #number of observations in test data
               x.train,   #p*n training data x
               y.train,   #n*1 training data y
               x.test,    #p*np test data x
               ntree,
               numcut,
               ndpost,
               nskip,
               power,
               base,
               tau,
               oldgamma,
               oldz,
               xinfo,
               offset,
               printevery)
  

    y_num <- length(unique(y.train))
    
    res$yhat.train  <- res$yhat.train  + offset
    #res$prob.train <- lapply(1:y_num, function(i){
    #  if(i==1){
    #    ## probability for outcome category 1
    #    apply(-res$yhat.train,2,pnorm)
    #  }
    #  else if(i==y_num){
    #    1 - apply(res$gammahat[,i] - res$yhat.train,2,pnorm)
    #  }
    #  else{
    #    apply(res$gammahat[,i+1] - res$yhat.train,2,pnorm) - apply(res$gammahat[,i] - res$yhat.train,2,pnorm)
    #  }
    #})
    
    res$offset <- offset
    
  if(np>0) {
    res$yhat.test <- res$yhat.test + offset
    #res$prob.test <- lapply(1:y_num, function(i){
    #  if(i==1){
    #    ## probability for outcome category 1
    #    apply(-res$yhat.test,2,pnorm)
    #  }
    #  else if(i==y_num){
    #    1 - apply(res$gammahat[,i] - res$yhat.test,2,pnorm)
    #  }
    #  else{
    #    apply(res$gammahat[,i+1] - res$yhat.test,2,pnorm) - apply(res$gammahat[,i] - res$yhat.test,2,pnorm)
    #  }
    #})
  } else {
    res$yhat.test <- NULL
    #res$prob.test <- NULL
  }
  
  attr(res, 'class') <- 'obart'
  return(res)
}
