
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

print.maxim <- function(x, hessian=FALSE) {
   cat("--------------------------------------------\n")
   cat(maximType(x), "\n")
   cat("Number of iterations:", nIter(x), "\n")
   cat("Return code:", returnCode(x), "\n")
   cat(returnMessage(x), "\n")
   if(length(x@lastStep@f0) != 0) {
      show(x@lastStep)
   }
   if(!is.null(x@estimate)) {
      cat("Function value:", x@maximum, "\n")
      out <- cbind(Estimate=x@estimate, Gradient=x@gradient, Active=activePar(x))
      print(out)
      if(hessian) {
         cat("Hessian:\n")
         print(x@hessian)
      }
   }
   cat("--------------------------------------------\n")
}
setMethod("print", "maxim", print.maxim)
setMethod("show", "maxim", function(object) print.maxim(object))

### (redundant) summary method for maxim.  Nothing to summarize, just print it
summary.maxim <- function(object) object
setMethod("summary", "maxim", summary.maxim)
rm(summary.maxim)

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
