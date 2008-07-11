## Return #of iterations for maxim objects

nIter.maxim <- function(x)
    x@iterations
setMethod("nIter", "maxim", nIter.maxim)
rm(nIter.maxim)
