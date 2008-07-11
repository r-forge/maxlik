library(methods)

setClass("lastStep",
         representation(theta0 = "numeric",
                        f0 = "numeric",
                        climb = "numeric"))

lastStep <- function(theta0, f0, climb) {
   new("lastStep", theta0, f0, climb)
}

showLastStep <- function(x) {
   cat("Function value:", x@f0)
   print(rbind(theta=x@theta0, climb=x@climb))
}
setMethod("show", "lastStep", showLastStep)
rm(showLastStep)

validLastStep <- function(object) {
   if(length(object@climb) != length(object@theta0))
       return(paste("'climb' and 'theta0' must be of equal length (currently ",
                    length(object@climb), " and ", length(object@theta0), ")"))
   if(length(object@f0) != 1)
       return(paste("'f0' must be scalar function value (currently ", lengh(object@f0), ")"))
   return(TRUE)
}
setValidity("lastStep", validLastStep)
rm(validLastStep)
