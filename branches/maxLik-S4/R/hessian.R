## Return Hessian of an object

hessian.maxim <- function(x)
    x@hessian
setMethod("hessian", "maxim", hessian.maxim)
rm(hessian.maxim)
