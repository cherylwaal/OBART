#' predict.obart
#' 
#'
#' @param object Returned from previous BART fit.
#' @param newdata Matrix of covariates to predict $y$.
#' @param mc.cores Number of threads to utilize.
#' @param openmp Logical value dictating whether OpenMP is utilized for parallel processing. Of course, this depends on whether OpenMP is available on your system which, by default, is verified with mc.cores.openmp.
#' @return What the function returns.
#' @export


predict.obart <- function(object, newdata, mc.cores=1, openmp=(mc.cores.openmp()>0), ...) {

    ##if(class(newdata) != "matrix") stop("newdata must be a matrix")

    #p <- length(object$treedraws$cutpoints)

    #if(p!=ncol(newdata))
    #    stop(paste0('The number of columns in newdata must be equal to ', p))

    if(.Platform$OS.type == "unix") mc.cores.detected <- detectCores()
    else mc.cores.detected <- NA

    call <- pwbart

    ##return(call(newdata, object$treedraws, mc.cores=mc.cores, mu=object$binaryOffset, ...))

    pred <- list(yhat.test=call(newdata, object$treedraws, mc.cores=mc.cores,
                                mu=object$offset, ...))
    
    y_num <- ncol(object$gammahat) - 1
    
    #pred$prob.test <- lapply(1:y_num, function(i){
    #  if(i==1){
    #    ## probability for outcome category 1
    #    apply(-pred$yhat.test,2,pnorm)
    #  }
    #  else if(i==y_num){
    #    1 - apply(object$gammahat[,i] - pred$yhat.test,2,pnorm)
    #  }
    #  else{
    #    apply(object$gammahat[,i+1] - pred$yhat.test,2,pnorm) - apply(object$gammahat[,i] - pred$yhat.test,2,pnorm)
    #  }
    #})
    
    pred$offset <- object$offset
    attr(pred, 'class') <- 'obart'

    return(pred)
}

