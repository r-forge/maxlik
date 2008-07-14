coef.maxim <- function( object) {
   return( object@estimate )
}
setMethod("coef", "maxim", coef.maxim)
rm(coef.maxim)
