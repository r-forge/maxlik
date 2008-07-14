### Methods for accessing loglik value maximum likelihood estimates

logLik.maxLik <- function(object)
    maxValue(object)
setMethod("logLik", "maxLik", logLik.maxLik)
rm(logLik.maxLik)
