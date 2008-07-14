## Return the #of parameters of model

nParam.lm <- function(x, ...)
    length(coefficients(x))

nParam.maxim <- function(x, free=FALSE, ...) {
   if(free)
       sum(activePar(x))
   else
       length(coef(x))
}
setMethod("nParam", "maxim", nParam.maxim)
rm(nParam.maxim)
