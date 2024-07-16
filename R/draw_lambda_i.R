

draw_lambda_i=function(lambda, mean, kmax=1000, thin=1)
    .Call("cdraw_lambda_i", lambda, mean, kmax, thin)

