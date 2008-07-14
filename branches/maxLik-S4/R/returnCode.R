## Returns return (error) code

returnCode.maxim <- function(x)
    x@code
setMethod("returnCode", "maxim", returnCode.maxim)
rm(returnCode.maxim)
