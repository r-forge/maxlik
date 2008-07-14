
returnMessage.maxim <- function(x)
    x@message
setMethod("returnMessage", "maxim", returnMessage.maxim)
rm(returnMessage.maxim)
