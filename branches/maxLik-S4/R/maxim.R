library(methods)
setClass("lastStep",
         representation(theta0 = "numeric",
                        f0 = "numeric",
                        climb = "numeric"))

setClass("maxim",
         representation(maximum = "numeric",
                        estimate = "numeric",
                        gradient = "vector",
                        hessian = "matrix",
                        code = "integer",
                        message = "character",
                        iterations = "integer",
                        lastStep = "lastStep",
                        activePar = "logical",
                        type = "character"))

maxim <- function(maximum, estimate, gradient, hessian=NULL,
                  code, message, iterations,
                  lastStep=NULL,
                  activePar=rep(TRUE, length(estimate)),
                  type) {
   new("maxim",
       maximum=maximum, estimate=estimate, gradient=gradient, hessian=hessian,
       code=code, message=message, iterations=iterations,
       lastStep=lastStep,
       activePar=activePar,
       type=type) 
}

lastStep <- function(theta0, f0, climb) {
   new("lastStep", theta0, f0, climb)
}

validMaxim <- function(object) {
   parLength <- length(object@estimate)
   if(length(object@activePar) != parLength)
       return(paste("'activePar' must have the same length as 'estimate' (currently",
                    length(object@activePar), "and", parLength, ")"))
   if(!is.null(object@hessian) &
      (nrow(object@hessian) != parLength | ncol(object@hessian) != parLength))
       return(paste("'hessian' must be either NULL or ", parLength, "x", parLength,
                    " matrix (currently ", nrow(object@hessian), "x", ncol(object@hessian), ")",
                    sep=""))
   return(TRUE)
}
setValidity("maxim", validMaxim)
rm(validMaxim)

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
