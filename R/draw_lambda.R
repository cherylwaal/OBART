

draw_lambda=function(lambda, mean, kmax=1000, thin=1)
    .Call("cdraw_lambda", lambda, mean, kmax, thin)

