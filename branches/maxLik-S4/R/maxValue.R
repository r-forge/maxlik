maxValue.maxim <- function(x)
    x@maximum
setMethod("maxValue", "maxim", maxValue.maxim)
rm(maxValue.maxim)
