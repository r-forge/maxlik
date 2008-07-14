
print.summary.maxLik <- function(x) {
   cat("--------------------------------------------\n")
   cat("Maximum Likelihood estimation\n")
   cat(maximType(x), ", ", nIter(x), " iterations\n", sep="")
   cat("Return code ", returnCode(x), ": ", returnMessage(x), "\n", sep="")
   if(!is.null(coef(x))) {
      cat("Log-Likelihood:", maxValue(x), "\n")
      cat(x@NActivePar, " free parameters\n")
      cat("Estimates:\n")
      printCoefmat(x@results)
   }
   cat("--------------------------------------------\n")
}
setMethod("show", "summary.maxLik", function(object) print.summary.maxLik(object))
setMethod("print", "summary.maxLik", print.summary.maxLik)

summary.maxLik <- function(object, results, NActivePar) {
         new("summary.maxLik",
             maximum= object@maximum,
             estimate= object@estimate,
             gradient= object@gradient,
             hessian= object@hessian,
             code= object@code,
             message= object@message,
             iterations= object@iterations,
             lastStep= object@lastStep,
             activePar= object@activePar,
             type= object@type,
             results=results,
             NActivePar=NActivePar)
}       

summaryMaxLik <- function( object, eigentol=1e-9,... ) {
   ## object      object of class "maxLik"
   ## 
   ## RESULTS:
   ## list of class "summary.maxLik" with following components:
   ## maximum    : function value at optimum
   ## estimate   : estimated parameter values at optimum
   ## gradient   :           gradient at optimum
   ## code       : code of convergence
   ## message    : message, description of the code
   ## iterations : number of iterations
   ## type       : type of optimisation
   ##
   nParam <- length(coef <- coef(object))
   activePar <- activePar(object)
   if(returnCode(object) < 100) {
      if(min(abs(eigen(hessian(object)[activePar,activePar],
                       symmetric=TRUE, only.values=TRUE)$values)) > eigentol) {
         varcovar <- matrix(0, nParam, nParam)
         varcovar[activePar,activePar] <-
             solve(-hessian(object)[activePar,activePar])
         hdiag <- diag(varcovar)
         if(any(hdiag < 0)) {
            warning("Diagonal of variance-covariance matrix not positive!\n")
         }
         stdd <- sqrt(hdiag)
         t <- coef/stdd
         p <- 2*pnorm( -abs( t))
      } else {
         stdd <- 0
         t <- 0
         p <- 0
      }
      results <- cbind("Estimate"=coef, "Std. error"=stdd, "t value"=t, "Pr(> t)"=p)
      Hess <- NULL
   } else {
      results <- NULL
      Hess <- NULL
   }
   summary <- summary.maxLik(object, results, sum(activePar))
   summary
}
setMethod("summary", "maxLik", summaryMaxLik)
rm(summaryMaxLik)
