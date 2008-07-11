maximTypeMaxim <- function(x)
    x@type
setMethod("maximType", "maxim", maximTypeMaxim)
rm(maximTypeMaxim)
