## Akaike (and other) information criteria
AIC.maxLik <- function(object, ..., k = 2)
    -2*logLik(object) + k*nParam(object, free=TRUE)
setMethod("AIC", "maxLik", AIC.maxLik)
rm(AIC.maxLik)
