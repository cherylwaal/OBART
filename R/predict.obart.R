#' predict.obart
#' 
#'
#' 
#' @param x.train Description of the first parameter.
#' @param y.train Description of the second parameter.
#' @param x.test Description of the second parameter.
#' @param xinfo Description of the second parameter.
#' @param usequants Description of the second parameter.
#' @param cont Description of the second parameter.
#' @param rm.const Description of the second parameter.
#' @param k Description of the second parameter.
#' @param power Description of the second parameter.
#' @param base Description of the second parameter.
#' @param rho Description of the second parameter.
#' @param ntree Description of the second parameter.
#' @param numcut Description of the second parameter.
#' @param ndpost Description of the second parameter.
#' @param nskip Description of the second parameter.
#' @param oldgamma Description of the second parameter.
#' @param oldz Description of the second parameter.
#' @param mc.cores Description of the second parameter.
#' @param nice Description of the second parameter.
#' @param seed Description of the second parameter.
#' 
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

