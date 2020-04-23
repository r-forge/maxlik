maxSGACompute <- function(fn,
                         start, 
                           # maximum lambda for Marquardt (1963)
                         finalHessian=FALSE,
                         bhhhHessian = FALSE,
                         fixed=NULL,
                         control=maxControl(),
                         ...) {
   ## Stochastic Gradient Ascent
   ## Parameters:
   ## fn          - the function to be maximized.  Returns either scalar or
   ##               vector value with possible attributes 
   ##               constPar and newVal
   ##               fn must return the value with attributes 'gradient'
   ##               and 'hessian'
   ##               fn must have an argument sumObs
   ## start       - initial parameter vector (eventually w/names)
   ## control       MaxControl object:
   ##     steptol     - minimum step size
   ##     lambda0       initial Hessian corrector (see Marquardt, 1963, p 438)
   ##     lambdaStep    how much Hessian corrector lambda is changed between
   ##                   two lambda trials
   ##                  (nu in Marquardt (1963, p 438)
   ##     maxLambda     largest possible lambda (if exceeded will give step error)
   ##     lambdatol   - max lowest eigenvalue when forcing pos. definite H
   ##     qrtol       - tolerance for qr decomposition
   ##     qac           How to handle the case where new function value is
   ##               smaller than the original one:
   ##                  "stephalving"   smaller step in the same direction
   ##                  "marquardt"     Marquardt (1963) approach
   ##     The stopping criteria
   ##     tol         - maximum allowed absolute difference between sequential values
   ##     reltol      - maximum allowed reltive difference (stops if < reltol*(abs(fn) + reltol)
   ##     gradtol     - maximum allowed norm of gradient vector
   ## 
   ##     iterlim     - maximum # of iterations
   ##     
   ## finalHessian  include final Hessian?  As computing final hessian does not carry any extra penalty for NR method, this option is
   ##               mostly for compatibility reasons with other maxXXX functions.
   ##               TRUE/something else  include
   ##               FALSE                do not include
   ## fixed       - a logical vector -- which parameters are taken as fixed.
   ##               Other paramters are treated as variable (free).
   ## ...           additional argument to 'fn'.  This may include
   ##               'fnOrig', 'gradOrig', 'hessOrig' if called fromm
   ##               'maxNR'.
   ##
   ## RESULTS:
   ## a list of class "maxim":
   ## maximum     function value at maximum
   ## estimate    the parameter value at maximum
   ## gradient        gradient
   ## hessian         Hessian
   ## code        integer code of success:
   ##             1 - gradient close to zero
   ##             2 - successive values within tolerance limit
   ##             3 - could not find a higher point (step error)
   ##             4 - iteration limit exceeded
   ##             5 - infinite function value
   ##             6  infinite gradient
   ##             7  infinite Hessian
   ##             100 - initial value out of range
   ## message     character message describing the code
   ## last.step   only present if code == 3 (step error).  A list with following components:
   ##             theta0    - parameter value which led to the error
   ##             f0        - function value at these parameter values
   ##             climb     - the difference between theta0 and the new approximated parameter value (theta1)
   ##             fixed     - logical vector, which parameters are constant (fixed, inactive, non-free)
   ## fixed       logical vector, which parameters were treated as constant (fixed, inactive, non-free)
   ## iterations  number of iterations
   ## type        "Newton-Raphson maximisation"
   ##
   ## References:
   ## Marquardt (1963), "An algorithm for least-squares estimation of nonlinear
   ##      parameters", J. Soc. Indust. Appl. Math 11(2), 431-441
   ##      
   ## -------------------------------------------------
   maximType <- "Stochastic Gradient Ascent"
   nimed <- names(start)
   start1 <- start
   nParam <- length(start1)
   learningRate <- slot(control, "SGA_learningRate")
   iter <- 0
   returnHessian <- ifelse( bhhhHessian, "BHHH", TRUE )
   f1 <- fn(start1, fixed = fixed, sumObs = TRUE, ...)
                           # have to compute fn as we cannot get gradient otherwise
   if(slot(control, "printLevel") > 0) {
      cat("Initial function value:", f1, "\n")
      if( isTRUE( attr( f1, "gradBoth" ) ) ) {
         warning( "the gradient is provided both as attribute 'gradient' and",
                 " as argument 'grad': ignoring argument 'grad'" )
      }
      if( isTRUE( attr( f1, "hessBoth" ) ) ) {
         warning( "the Hessian is provided both as attribute 'hessian' and",
                 " as argument 'hess': ignoring argument 'hess'" )
      }
   }
   G1 <- attr( f1, "gradient" )
   if(slot(control, "printLevel") > 2) {
      cat("Initial gradient value:\n")
      print(G1)
   }
   if(any(is.na(G1[!fixed]))) {
      stop("NA in the initial gradient")
   }
   if(any(is.infinite(G1[!fixed]))) {
      stop("Infinite initial gradient")
   }
   if(length(G1) != nParam) {
      stop( "length of gradient (", length(G1),
         ") not equal to the no. of parameters (", nParam, ")" )
   }
   if( slot(control, "printLevel") > 1) {
      cat( "----- Initial parameters: -----\n")
      cat( "fcn value:",
      as.vector(f1), "\n")
      a <- cbind(start1, G1, as.integer(!fixed))
      dimnames(a) <- list(nimed, c("parameter", "initial gradient",
                                          "free"))
      print(a)
   }
   ## ---------------- Main interation loop ------------------------
   ## we do not need to compute the function itself here, except for
   ## printing
   repeat {
      ## break here if iterlim == 0
      if( iter >= slot(control, "iterlim")) {
         code <- 4; break
      }
      iter <- iter + 1
      start0 <- start1
      f0 <- f1
      G0 <- G1
      if(any(is.na(G0[!fixed]))) {
         stop("NA in gradient")
      }
      start1 <- start0 + learningRate*G0
      if( slot(control, "printLevel") > 1) {
        cat( "-----Iteration", iter, "-----\n")
      }
      if( iter >= slot(control, "iterlim")) {
         code <- 4; break
      }
      ## break here to avoid potentially costly gradient computation
      if( iter >= slot(control, "iterlim")) {
         code <- 4; break
      }
      ## still iterations to go, hence compute gradient
      f1 <- fn(start1, fixed = fixed, sumObs = TRUE,
               returnHessian = returnHessian, ...)
                           # The call calculates new function,
                           # gradient, and Hessian values
      G1 <- attr( f1, "gradient" )
      if(any(is.na(G1[!fixed]))) {
         cat("Iteration", iter, "\n")
         cat("Parameter:\n")
         print(start1)
         print(head(G1, n=30))
         stop("NA in gradient")
      }
      if(any(is.infinite(G1))) {
         code <- 6; break;
      }
      if(slot(control, "printLevel") > 2) {
         cat(" learning rate", learningRate, " fcn value:",
            formatC(as.vector(f1), digits=8, format="f"),  "\n")
         a <- cbind(learningRate*G0, start1, G1, as.integer(!fixed))
         dimnames(a) <- list(names(start0), c("amount", "new param",
                                             "new gradient", "active"))
         print(a)
      }
      if( sqrt( crossprod( G1[!fixed] ) ) < slot(control, "gradtol") ) {
         code <-1; break
      }
   }  # main iteration loop
   if( slot(control, "printLevel") > 0) {
      cat( "--------------\n")
      cat( maximMessage( code), "\n")
      cat( iter, " iterations\n")
      cat( "estimate:", start1, "\n")
      cat( "Function value:", f1, "\n")
   }
   names(start1) <- nimed
   if(finalHessian & !bhhhHessian) {
      F1 <- fn( start1, fixed = fixed, sumObs = FALSE,
               returnHessian = TRUE, ... )
      G1 <- attr( F1, "gradient" )
   }
   if(observationGradient(G1, length(start1))) {
      gradientObs <- G1
      colnames( gradientObs ) <- nimed 
      G1 <- colSums(as.matrix(G1 ))
   }
   else {
      gradientObs <- NULL
   }
   names( G1 ) <- nimed
   ## calculate (final) Hessian
   if(tolower(finalHessian) == "bhhh") {
      if(!is.null(gradientObs)) {
         hessian <- - crossprod( gradientObs )
         attr(hessian, "type") <- "BHHH"
      } else {
         hessian <- NULL
         warning("For computing the final Hessian by 'BHHH' method, the log-likelihood or gradient must be supplied by observations")
      }
   } else if( finalHessian != FALSE ) {
      hessian <- attr( F1, "hessian" )
   } else {
       hessian <- NULL
   }
   if( !is.null( hessian ) ) {
      rownames( hessian ) <- colnames( hessian ) <- nimed
   }

   ## remove attributes from final value of objective (likelihood) function
   attributes( f1 )$gradient <- NULL
   attributes( f1 )$hessian <- NULL
   attributes( f1 )$gradBoth <- NULL
   attributes( f1 )$hessBoth <- NULL
   ##
   result <- list(
       maximum = unname( drop( f1 ) ),
                  estimate=start1,
                  gradient=drop(G1),
      hessian=hessian,
                  code=code,
       message=maximMessage( code),
      fixed=fixed,
       iterations=iter,
                  type=maximType)
   if( exists( "gradientObs" ) ) {
      result$gradientObs <- gradientObs
   }
   result <- c(result, control=control)
                           # attach the control parameters
   ##
   class(result) <- c("maxim", class(result))
   invisible(result)
}

returnCode.maxim <- function(x, ...)
    x$code
