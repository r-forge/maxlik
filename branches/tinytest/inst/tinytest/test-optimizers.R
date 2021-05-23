### This code tests all the methods and main parameters.  It includes:
### * analytic gradients/Hessian
### * fixed parameters
### * inequality constraints
### * equality constraints

## do not run unless 'NOT_CRAN' explicitly defined
## (Suggested by Sebastian Meyer and others)
if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
    message("skipping slow optimizer tests")
    q("no")
}
if(!requireNamespace("tinytest", quietly = TRUE)) {
   message("These tests require 'tinytest' package\n")
   q("no")
}
library(maxLik)

## data to fit a normal distribution
# set seed for pseudo random numbers
set.seed( 123 )
tol <- .Machine$double.eps^0.25
## generate a variable from normally distributed random numbers
truePar <- c(mu=1, sigma=2)
x <- rnorm( 100, truePar[1], truePar[2] )
xSaved <- x

## log likelihood function
llf <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   sum(dnorm(x, mu, sigma, log=TRUE))
}

## log likelihood function (individual observations)
llfInd <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   llValues <- -0.5 * log( 2 * pi ) - log( sigma ) -
      0.5 * ( x - mu )^2 / sigma^2
   return( llValues )
}

## function to calculate analytical gradients
gf <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   N <- length( x )
   llGrad <- c( sum( ( x - mu ) / sigma^2 ),
      - N / sigma + sum( ( x - mu )^2 / sigma^3 ) )
   return( llGrad )
}

## function to calculate analytical gradients (individual observations)
gfInd <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   llGrads <- cbind( ( x - mu ) / sigma^2,
      - 1 / sigma + ( x - mu )^2 / sigma^3 )
   return( llGrads )
}

## log likelihood function with gradients as attributes
llfGrad <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   N <- length( x )
   llValue <- -0.5 * N * log( 2 * pi ) - N * log( sigma ) -
      0.5 * sum( ( x - mu )^2 / sigma^2 )
   attributes( llValue )$gradient <- c( sum( ( x - mu ) / sigma^2 ),
      - N / sigma + sum( ( x - mu )^2 / sigma^3 ) )
   return( llValue )
}

## log likelihood function with gradients as attributes (individual observations)
llfGradInd <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   llValues <- -0.5 * log( 2 * pi ) - log( sigma ) -
      0.5 * ( x - mu )^2 / sigma^2
   attributes( llValues )$gradient <- cbind( ( x - mu ) / sigma^2,
      - 1 / sigma + ( x - mu )^2 / sigma^3 )
   return( llValues )
}

## function to calculate analytical Hessians
hf <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   N <- length( x )
   llHess <- matrix( c(
      N * ( - 1 / sigma^2 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      N / sigma^2 + sum( - 3 * ( x - mu )^2 / sigma^4 ) ),
      nrow = 2, ncol = 2 )
   return( llHess )
}

## log likelihood function with gradients and Hessian as attributes
llfGradHess <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   N <- length( x )
   llValue <- -0.5 * N * log( 2 * pi ) - N * log( sigma ) -
      0.5 * sum( ( x - mu )^2 / sigma^2 )
   attributes( llValue )$gradient <- c( sum( ( x - mu ) / sigma^2 ),
      - N / sigma + sum( ( x - mu )^2 / sigma^3 ) )
   attributes( llValue )$hessian <- matrix( c(
      N * ( - 1 / sigma^2 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      N / sigma^2 + sum( - 3 * ( x - mu )^2 / sigma^4 ) ),
      nrow = 2, ncol = 2 )
   return( llValue )
}

## log likelihood function with gradients as attributes (individual observations)
llfGradHessInd <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   N <- length( x )
   llValues <- -0.5 * log( 2 * pi ) - log( sigma ) -
      0.5 * ( x - mu )^2 / sigma^2
   attributes( llValues )$gradient <- cbind( ( x - mu ) / sigma^2,
      - 1 / sigma + ( x - mu )^2 / sigma^3 )
   attributes( llValues )$hessian <- matrix( c(
      N * ( - 1 / sigma^2 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      N / sigma^2 + sum( - 3 * ( x - mu )^2 / sigma^4 ) ),
      nrow = 2, ncol = 2 )
   return( llValues )
}


# start values
startVal <- c( mu = 0, sigma = 1 )

## basic NR: test if all methods work
ml <- maxLik( llf, start = startVal )
expect_equal(
   coef(ml), truePar, tol=2*max(stdEr(ml))
)
expect_stdout(
   print( ml ),
   pattern = "Estimate\\(s\\): 1.18.*1.81"
)
expect_stdout(
   print( summary( ml )),
   pattern = "Estimates:"
)
expect_equal(
   activePar( ml ), c(mu=TRUE, sigma=TRUE)
)
expect_equal(
   AIC( ml ), 407.167892384587,
   tol = 0.1, check.attributes=FALSE
)
expect_equal(
   coef( ml ), c(mu=1.181, sigma=1.816),
   tol = 0.001
)
expect_stdout(
   condiNumber( ml, digits = 3),
   "mu[[:space:]]+1[[:space:]\n]+sigma[[:space:]]+1\\."
)
expect_equal(
   hessian( ml), matrix(c(-30.3, 0, 0, -60.6), 2, 2),
   tol = 0.01, check.attributes = FALSE
)
expect_equal(
   logLik( ml ), -201.583946192294,
   tol = tol, check.attributes = FALSE
)
expect_equal(
   maximType( ml ), "Newton-Raphson maximisation"
)
expect_equal(
   nIter( ml ) > 5, TRUE
)
expect_error(
   nObs( ml ),
   "cannot return the number of observations"
)
expect_equal(
   nParam( ml ), 2
)
expect_equal(
   returnCode( ml ), 1
)
expect_equal(
   returnMessage( ml ), "gradient close to zero (gradtol)"
)
expect_equal(
   vcov( ml ), matrix(c(0.032975, 0, 0, 0.0165), 2, 2),
   tol=0.01, check.attributes = FALSE
)
expect_equal(
   logLik( summary( ml ) ), logLik(ml)
)
mlInd <- maxLik( llfInd, start = startVal )
expect_stdout(
   print( summary( mlInd ), digits = 2 ),
   "mu +1\\.18"
)
expect_equal(
   nObs( mlInd ), length(x)
)
## Marquardt (1963) correction
mlM <- maxLik( llf, start = startVal, qac="marquardt")
expect_equal(
   coef(mlM), coef(ml),
                           # coefficients should be the same as above
   tol=tol
)
expect_equal(
   returnMessage(mlM), returnMessage(ml)
)

## with analytical gradients
## compare coefficients, Hessian
mlg <- maxLik( llf, gf, start = startVal )
expect_equal(coef(ml), coef(mlg), tol=tol)
expect_equal(hessian(ml), hessian(mlg), tolerance = 1e-3)
## gradient with individual components
mlgInd <- maxLik( llfInd, gfInd, start = startVal )
expect_equal(coef(mlInd), coef(mlgInd), tolerance = 1e-3)
expect_equal(hessian(mlg), hessian(mlgInd), tolerance = 1e-3)

## with analytical gradients as attribute
mlG <- maxLik( llfGrad, start = startVal )
expect_equal(coef(mlG), coef(mlg), tolerance = tol)
expect_equivalent(gradient(mlG), gf( coef( mlG ) ), tolerance = tol)
mlGInd <- maxLik( llfGradInd, start = startVal )
expect_equal(coef(mlGInd), coef(mlgInd), tolerance = tol)
expect_equivalent(gradient(mlGInd), colSums( gfInd( coef( mlGInd ) ) ), tolerance = tol)
expect_equivalent(estfun(mlGInd), gfInd( coef( mlGInd ) ), tolerance=tol)

## with analytical gradients as argument and attribute
expect_warning(mlgG <- maxLik( llfGrad, gf, start = startVal))
expect_equal(coef(mlgG), coef(mlg), tolerance = tol)

## with analytical gradients and Hessians
mlgh <- maxLik( llf, gf, hf, start = startVal )
expect_equal(coef(mlg), coef(mlgh), tolerance = tol)

## with analytical gradients and Hessian as attribute
mlGH <- maxLik( llfGradHess, start = startVal )
expect_equal(coef(mlGH), coef(mlgh), tolerance = tol)

## with analytical gradients and Hessian as argument and attribute
expect_warning(mlgGhH <- maxLik( llfGradHess, gf, hf, start = startVal ))
expect_equal(coef(mlgGhH), coef(mlgh), tolerance = tol)


## ---------- BHHH method ----------
## cannot do BHHH if llf not provided by individual
x <- xSaved[1]
expect_error( maxLik( llfInd, start = startVal, method = "BHHH" ) )
## 2 observations: can do BHHH
x <- xSaved[1:2]
expect_silent( maxLik( llfInd, start = startVal, method = "BHHH" ) )
##
x <- xSaved
mlBHHH <- maxLik( llfInd, start = startVal, method = "BHHH" )
expect_stdout(print( mlBHHH ),
              pattern = "Estimate\\(s\\): 1\\.18.* 1\\.81")
expect_stdout(print(summary( mlBHHH)), pattern = "mu *1.18")
expect_equivalent(activePar( mlBHHH ), c(TRUE, TRUE))
expect_equivalent(AIC( mlBHHH ), 407.168, tolerance=0.01)
expect_equal(coef( mlBHHH ), setNames(c(1.180808, 1.816485), c("mu", "sigma")), tolerance=tol)
expect_equal(condiNumber( mlBHHH, printLevel=0),
             setNames(c(1, 1.72), c("mu", "sigma")), tol=0.01)
expect_equivalent(hessian( mlBHHH ),
                  matrix(c(-30.306411, -1.833632, -1.833632, -55.731646), 2, 2),
                  tolerance=0.01)
expect_equivalent(logLik( mlBHHH ), -201.583946192983, tolerance=tol)
expect_equal(maximType( mlBHHH ), "BHHH maximisation")
expect_equal(nIter(mlBHHH) > 3, TRUE)
                           # here 12 iterations
expect_equal(nParam( mlBHHH ), 2)
expect_equal(returnCode( mlBHHH ), 8)
expect_equal(returnMessage( mlBHHH ),
             "successive function values within relative tolerance limit (reltol)")
expect_equivalent(vcov( mlBHHH ),
                  matrix(c(0.03306213, -0.00108778, -0.00108778, 0.01797892), 2, 2),
                  tol=0.001)
expect_equivalent(logLik(summary(mlBHHH)), -201.583946192983, tolerance=tol)
expect_equal(coef(ml), coef(mlBHHH), tol=tol)
expect_equal(stdEr(ml), stdEr(mlBHHH), tol=0.1)
expect_equal(nObs( mlBHHH ), length(x))
# final Hessian = usual Hessian
expect_silent(mlBhhhH <- maxLik( llfInd, start = startVal, method = "BHHH", 
                                finalHessian = TRUE )
              )
                           # do not test Hessian equality--BHHH may be imprecise, at least
                           # for diagonal elements
expect_stdout(print(hessian( mlBhhhH )),
              pattern="mu.*\nsigma.+")
## Marquardt (1963) correction
expect_silent(mlBHHHM <- maxLik( llfInd, start = startVal, method = "BHHH", qac="marquardt"))
expect_equal(coef(mlBHHHM), coef(mlBHHH), tolerance=tol)
expect_equal(returnMessage(mlBHHHM), "successive function values within relative tolerance limit (reltol)")

## BHHH with analytical gradients
expect_error( maxLik( llf, gf, start = startVal, method = "BHHH" ) )
                           # need individual log-likelihood
expect_error( maxLik( llfInd, gf, start = startVal, method = "BHHH" ) )
                           # need individual gradient
x <- xSaved[1]  # test with a single observation
expect_error(maxLik( llf, gfInd, start = startVal, method = "BHHH" ))
                           # gradient must have >= 2 rows
expect_error( maxLik( llfInd, gfInd, start = startVal, method = "BHHH" ) )
                           # ditto even if individual likelihood components
x <- xSaved[1:2]  # test with 2 observations
expect_silent(maxLik( llf, gfInd, start = startVal, method = "BHHH",
                     iterlim=1))
                           # should work with 2 obs
expect_silent( maxLik( llfInd, gfInd, start = startVal, method = "BHHH",
                      iterlim=1) )
                           # should work with 2 obs
x <- xSaved
expect_silent(mlgBHHH <- maxLik( llfInd, gfInd, start = startVal, method = "BHHH" ))
                           # individual log-likelihood, gradient
expect_equal(coef(mlBHHH), coef(mlgBHHH), tolerance = tol)
expect_equal(coef(mlg), coef(mlgBHHH), tolerance = tol)
expect_silent(mlgBHHH2 <- maxLik( llf, gfInd, start = startVal, method = "BHHH" ))
                           # aggregated log-likelihood, individual gradient
expect_equal(coef(mlgBHHH), coef(mlgBHHH2), tolerance=tol)
                           # final Hessian = usual Hessian
expect_silent(
   mlgBhhhH <- maxLik( llf, gfInd, start = startVal, method = "BHHH", 
                      finalHessian = TRUE )
)
expect_equal(hessian(mlgBhhhH), hessian(mlBhhhH), tolerance = 1e-2)

## with analytical gradients as attribute
expect_error( maxLik( llfGrad, start = startVal, method = "BHHH" ) )
                           # no individual gradients provided
x <- xSaved[1]
expect_error( maxLik( llfGrad, start = startVal, method = "BHHH" ),
             pattern = "gradient is not a matrix")
                           # get an error about need a matrix
expect_error( maxLik( llfGradInd, start = startVal, method = "BHHH" ),
             pattern = "at least as many rows")
                           # need at least two obs
x <- xSaved[1:2]
expect_error( maxLik( llfGrad, start = startVal, method = "BHHH" ),
             pattern = "gradient is not a matrix")
                           # enough obs but no individual grad
x <- xSaved
expect_silent(mlGBHHH <- maxLik( llfGradInd, start = startVal, method = "BHHH" ))
expect_equal(coef(mlGBHHH), coef(mlgBHHH), tolerance = tol)
                           # final Hessian = usual Hessian
expect_silent(mlGBhhhH <- maxLik( llfGradInd, start = startVal, method = "BHHH", 
                                 finalHessian = TRUE ))
expect_equal(hessian(mlGBhhhH), hessian(mlgBhhhH), tolerance = tol)

## with analytical gradients as argument and attribute
expect_warning(mlgGBHHH <- maxLik( llfGradInd, gfInd, start = startVal, method = "BHHH" ),
               pattern = "both as attribute 'gradient' and as argument 'grad'")
                           # warn about double gradient
expect_equal(coef(mlgGBHHH), coef(mlgBHHH), tolerance = tol)
## with unused Hessian
expect_silent(mlghBHHH <- maxLik( llfInd, gfInd, hf, start = startVal, method = "BHHH" ))
expect_equal(coef(mlgBHHH), coef(mlghBHHH), tolerance = tol)
## final Hessian = usual Hessian
expect_silent(
   mlghBhhhH <- maxLik( llfInd, gfInd, hf, start = startVal, method = "BHHH", 
                       finalHessian = TRUE )
)
expect_equivalent(hessian(mlghBhhhH), hessian(mlghBHHH), tolerance = 0.2)
                           # BHHH and ordinary hessian differ quite a bit
## with unused Hessian as attribute
expect_silent(mlGHBHHH <- maxLik( llfGradHessInd, start = startVal, method = "BHHH" ))
expect_equal(coef(mlGHBHHH), coef(mlghBHHH), tolerance = tol)
## final Hessian = usual Hessian
expect_silent(mlGHBhhhH <- maxLik( llfGradHessInd, start = startVal, method = "BHHH", 
                                  finalHessian = TRUE ))
expect_equal(hessian(mlGHBhhhH), hessian(mlghBhhhH), tolerance = tol)
## with analytical gradients and Hessian as argument and attribute
expect_warning(
   mlgGhHBHHH <- maxLik( llfGradHessInd, gfInd, hf, start = startVal, method = "BHHH" ),
   pattern = "both as attribute 'gradient' and as argument 'grad': ignoring"
)
expect_equal(coef(mlgGhHBHHH), coef(mlghBHHH), tolerance = tol)
expect_equal(hessian(mlgGhHBHHH), hessian(mlGHBHHH), tolerance = tol)

## ---------- Test BFGS methods ----------
optimizerNames <- c(bfgsr = "BFGSR", bfgs = "BFGS", nm = "Nelder-Mead",
                    sann = "SANN", cg = "CG")
successCodes <- list(bfgsr = 1:4, bfgs = 0, nm = 0, sann = 0, cg = 0)
successMsgs <- list(bfgsr = c("successive function values within tolerance limit (tol)"),
                    bfgs = c("successful convergence "),
                           # includes space at end...
                    nm = c("successful convergence "),
                    sann = c("successful convergence "),
                    cg = c("successful convergence ")
                    )
for(optimizer in c("bfgsr", "bfgs", "nm", "sann", "cg")) {
   expect_silent(mlResult <- maxLik( llf, start = startVal, method = optimizer ))
   expect_stdout(print( mlResult ),
                 pattern = paste0(optimizerNames[optimizer], " maximization")
                 )
   expect_stdout(print( summary( mlResult )),
                 pattern = paste0(optimizerNames[optimizer], " maximization,.*Estimates:")
                 )
   expect_equal(coef(ml), coef(mlResult), tolerance=0.001)
   expect_equal(stdEr(ml), stdEr(mlResult), tolerance=0.01)
   expect_equal(activePar( mlResult ), c(mu=TRUE, sigma=TRUE))
   expect_equivalent(AIC( mlResult ), 407.167893392749, tolerance=tol)
   expect_equivalent( hessian( mlResult ),
                     matrix(c(-30.32596, 0.00000, 0.00000, -60.59508), 2, 2),
                     tolerance = 0.01)
   expect_equivalent(logLik( mlResult ), -201.5839, tolerance = 0.01)
   expect_equal(maximType( mlResult ),
                paste0(optimizerNames[optimizer], " maximization")
                )
   expect_true(nIter( mlResult ) > 1 & is.integer(nIter(mlResult)))
   expect_error( nObs( mlResult ),
                pattern = "cannot return the number of observations")
   expect_equal(nParam( mlResult ), 2)
   expect_true(returnCode( mlResult ) %in% successCodes[[optimizer]])
   expect_equal(returnMessage( mlResult), successMsgs[[optimizer]])
   expect_equal(logLik( summary( mlResult ) ), logLik(mlResult))
   ## individual observations
   expect_silent(mlIndResult <- maxLik( llfInd, start = startVal, method = optimizer))
   expect_stdout(print( summary( mlIndResult )),
                 pattern = paste0(optimizerNames[optimizer], " maximization,.*Estimates:")
                 )
   expect_equal(coef(mlResult), coef(mlIndResult), tolerance = tol)
   expect_equal(stdEr(mlResult), stdEr(mlIndResult), tolerance = 0.01)
   expect_equal(nObs( mlIndResult ), length(x))
   ## with analytic gradients
   expect_silent(mlgResult <- maxLik( llf, gf, start = startVal, method = optimizer))
   expect_equal(coef(mlgResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlgResult), stdEr(mlResult), tolerance = 0.01)
   expect_silent(mlgIndResult <- maxLik( llfInd, gfInd, start = startVal,
                                        method = optimizer ))
   expect_equal(coef(mlgIndResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlgIndResult), stdEr(mlResult), tolerance = 0.01)
   ## with analytical gradients as attribute
   expect_silent(mlGResult <- maxLik( llfGrad, start = startVal,
                                     method = optimizer))
   expect_equal(coef(mlGResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlGResult), stdEr(mlResult), tolerance = 0.01)
   expect_silent(mlGIndResult <- maxLik( llfGradInd, start = startVal, method = optimizer ))
   expect_equal(coef(mlGIndResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlGIndResult), stdEr(mlResult), tolerance = 0.01)
   ## with analytical gradients as argument and attribute
   expect_warning(mlgGResult <- maxLik( llfGrad, gf, start = startVal, method = optimizer ))
   expect_equal(coef(mlgGResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlgGResult), stdEr(mlResult), tolerance = 0.01)
   ## with analytical gradients and Hessians
   expect_silent(mlghResult <- maxLik( llf, gf, hf, start = startVal, method = optimizer ))
   expect_equal(coef(mlghResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlghResult), stdEr(mlResult), tolerance = 0.01)
   ## with analytical gradients and Hessian as attribute
   expect_silent(mlGHResult <- maxLik( llfGradHess, start = startVal, method = optimizer ))
   expect_equal(coef(mlGHResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlGHResult), stdEr(mlResult), tolerance = 0.01)
   ## with analytical gradients and Hessian as argument and attribute
   expect_warning(mlgGhHResult <- maxLik( llfGradHess, gf, hf, start = startVal, method = optimizer ))
   expect_equal(coef(mlgGhHResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlgGhHResult), stdEr(mlResult), tolerance = 0.01)
}


### ---------- with fixed parameters ----------
## start values
startValFix <- c( mu = 1, sigma = 1 )
## fix mu (the mean ) at its start value
isFixed <- c( TRUE, FALSE )
successMsgs <- list(bfgsr = c("successive function values within tolerance limit (tol)"),
                    bfgs = c("successful convergence "),
                           # includes space at end...
                    nm = c("successful convergence "),
                    sann = c("successful convergence "),
                    cg = c("successful convergence ")
                    )
## NR method with fixed parameters
for(optimizer in c("nr", "bfgsr", "bfgs", "nm", "sann", "cg")) {
   expect_silent(
      mlFix <- maxLik( llf, start = startValFix, fixed = isFixed, method=optimizer)
   )
   expect_equivalent(coef(mlFix)[1], 1)
   expect_equivalent(stdEr(mlFix)[1], 0)
   cat(optimizer, "\n")
   print(hessian(mlFix))
   expect_equivalent(hessian( mlFix),
                     matrix(c(NA, NA, NA, -59.99823), 2, 2),
                     tolerance=0.01)
   mlFix3 <- maxLik(llf, start = startValFix, fixed = "mu", method=optimizer)
   expect_equal(coef(mlFix), coef(mlFix3))
   mlFix4 <- maxLik( llf, start = startValFix, fixed = which(isFixed),
                    method=optimizer)
   expect_equal(coef(mlFix), coef(mlFix4))
   expect_equivalent(activePar( mlFix ), !isFixed)
   expect_equal(nParam( mlFix ), 2)
   ## with analytical gradients
   mlgFix <- maxLik( llf, gf, start = startValFix, fixed = isFixed,
                    method=optimizer)
   expect_equal(coef(mlgFix), coef(mlFix), tolerance=tol)
   ## with analytical gradients and Hessians
   mlghFix <- maxLik( llf, gf, hf, start = startValFix, fixed = isFixed,
                     method=optimizer)
   expect_equal(coef(mlghFix), coef(mlFix), tolerance=tol)
}

## BHHH method with fixed parameters--need individual obs
mlFixBHHH <- maxLik( llfInd, start = startValFix, fixed = isFixed,
                    method = "BHHH" )
expect_equal(coef(mlFixBHHH), coef(mlFix))
ar
## mlFixBHHH1 <- maxLik( llfInd, start = startValFix, activePar = 2,
##    method = "BHHH" )
## expect_equal( mlFixBHHH, mlFixBHHH1, tolerance = 1e-3 )
## mlFixBHHH2 <- maxLik( llfInd, start = startValFix, fixed = isFixed,
##    method = "BHHH" )
## expect_equal( mlFixBHHH, mlFixBHHH2, tolerance = 1e-3 )
## mlFixBHHH3 <- maxLik( llfInd, start = startValFix, fixed = "mu",
##    method = "BHHH" )
## expect_equal( mlFixBHHH, mlFixBHHH3, tolerance = 1e-3 )
## mlFixBHHH4 <- maxLik( llfInd, start = startValFix, fixed = 1,
##    method = "BHHH" )
## expect_equal( mlFixBHHH, mlFixBHHH4, tolerance = 1e-3 )
## print( mlFixBHHH )
## print( summary( mlFixBHHH ), digits = 2 )
## activePar( mlFixBHHH )
## AIC( mlFixBHHH )
## coef( mlFixBHHH )
## condiNumber( mlFixBHHH, digits = 3 )
## round( hessian( mlFixBHHH ), 1 )
## logLik( mlFixBHHH )
## maximType( mlFixBHHH )
## nIter( mlFixBHHH )
## nParam( mlFixBHHH )
## returnCode( mlFixBHHH )
## returnMessage( mlFixBHHH )
## round( vcov( mlFixBHHH ), 3 )
## logLik( summary( mlFixBHHH ) )
## expect_equal( mlFix[ -c( 4, 5, 6, 9, 10 ) ], mlFixBHHH[ -c( 4, 5, 6, 9, 10, 11 ) ],
##    tolerance = 1e-3 )
## round( mlFixBHHH[[ 11 ]], 3 )
## nObs( mlFixBHHH )

## # with analytical gradients
## mlgFixBHHH <- maxLik( llfInd, gfInd, start = startValFix, activePar = !isFixed,
##    method = "BHHH" )
## mlgFixBHHH1 <- maxLik( llfInd, gfInd, start = startValFix, activePar = 2,
##    method = "BHHH" )
## expect_equal( mlgFixBHHH, mlgFixBHHH1, tolerance = 1e-3 )
## mlgFixBHHH2 <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
##    method = "BHHH" )
## expect_equal( mlgFixBHHH, mlgFixBHHH2, tolerance = 1e-3 )
## mlgFixBHHH3 <- maxLik( llfInd, gfInd, start = startValFix, fixed = "mu",
##    method = "BHHH" )
## expect_equal( mlgFixBHHH, mlgFixBHHH3, tolerance = 1e-3 )
## mlgFixBHHH4 <- maxLik( llfInd, gfInd, start = startValFix, fixed = 1,
##    method = "BHHH" )
## expect_equal( mlgFixBHHH, mlgFixBHHH4, tolerance = 1e-3 )
## print( summary( mlgFixBHHH ), digits = 2 )
## expect_equal( mlFixBHHH, mlgFixBHHH, tolerance = 1e-3 )
## mlgFixBHHH2 <- maxLik( llf, gfInd, start = startValFix, activePar = !isFixed,
##    method = "BHHH")
## expect_equal( mlgFixBHHH, mlgFixBHHH2, tolerance = 1e-3 )

## # with unused Hessians
## mlghFixBHHH <- maxLik( llfInd, gfInd, hf, start = startValFix, activePar = !isFixed,
##    method = "BHHH" )
## expect_equal( mlgFixBHHH, mlghFixBHHH, tolerance = 1e-3 )

## ## BFGS method with fixed parameters
## mlFixBfgs <- maxLik( llf, start = startValFix, fixed = isFixed,
##    method = "BFGS" )
## mlFixBfgs3 <- maxLik( llf, start = startValFix, fixed = "mu",
##    method = "BFGS" )
## expect_equal( mlFixBfgs, mlFixBfgs3, tolerance = 1e-3 )
## mlFixBfgs4 <- maxLik( llf, start = startValFix, fixed = 1,
##    method = "BFGS" )
## expect_equal( mlFixBfgs, mlFixBfgs4, tolerance = 1e-3 )
## print( mlFixBfgs )
## print( summary( mlFixBfgs ), digits = 2 )
## activePar( mlFixBfgs )
## AIC( mlFixBfgs )
## coef( mlFixBfgs )
## condiNumber( mlFixBfgs, digits = 3 )
## round( hessian( mlFixBfgs ), 1 )
## logLik( mlFixBfgs )
## maximType( mlFixBfgs )
## nIter( mlFixBfgs )
## nParam( mlFixBfgs )
## returnCode( mlFixBfgs )
## returnMessage( mlFixBfgs )
## round( vcov( mlFixBfgs ), 3 )
## logLik( summary( mlFixBfgs ) )
## expect_equal( mlghFix[ -c( 5, 6, 9, 10 ) ], mlFixBfgs[ -c( 5, 6, 9, 10, 11 ) ],
##    tolerance = 1e-3 )
## mlIndFixBfgs <- maxLik( llfInd, start = startValFix, fixed = isFixed,
##    method = "BFGS" )
## expect_equal( mlFixBfgs[-c(4,9)], mlIndFixBfgs[ -c(4,9,12) ], tolerance = 1e-3 )
## print(formatC(mlIndFixBfgs$gradientObs, format="f", digits=4, width=7), quote=FALSE)
##                            # print fradient, only 4 digits to avoid clutter in R CMD tests
## mlIndFixBfgs3 <- maxLik( llfInd, start = startValFix, fixed = "mu",
##    method = "BFGS" )
## expect_equal( mlIndFixBfgs, mlIndFixBfgs3, tolerance = 1e-3 )
## mlIndFixBfgs4 <- maxLik( llfInd, start = startValFix, fixed = 1,
##    method = "BFGS" )
## expect_equal( mlIndFixBfgs, mlIndFixBfgs4, tolerance = 1e-3 )
## nObs( mlIndFixBfgs )

## # with analytical gradients
## mlgFixBfgs <- maxLik( llf, gf, start = startValFix, fixed = isFixed,
##    method = "BFGS" )
## mlgFixBfgs3 <- maxLik( llf, gf, start = startValFix, fixed = "mu",
##    method = "BFGS" )
## expect_equal( mlgFixBfgs, mlgFixBfgs3, tolerance = 1e-3 )
## mlgFixBfgs4 <- maxLik( llf, gf, start = startValFix, fixed = 1,
##    method = "BFGS" )
## expect_equal( mlgFixBfgs, mlgFixBfgs4, tolerance = 1e-3 )
## print( summary( mlgFixBfgs ), digits = 2 )
## expect_equal( mlFixBfgs[ -9 ], mlgFixBfgs[ -9 ], tolerance = 1e-3 )
## mlgIndFixBfgs <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
##    method = "BFGS")
## expect_equal( mlgFixBfgs[ ], mlgIndFixBfgs[ -12 ], tolerance = 1e-3 )
## round( mlgIndFixBfgs[[ 12 ]], 3 )
## mlgIndFixBfgs3 <- maxLik( llfInd, gfInd, start = startValFix, fixed = "mu",
##    method = "BFGS" )
## expect_equal( mlgIndFixBfgs, mlgIndFixBfgs3, tolerance = 1e-3 )
## mlgIndFixBfgs4 <- maxLik( llfInd, gfInd, start = startValFix, fixed = 1,
##    method = "BFGS" )
## expect_equal( mlgIndFixBfgs, mlgIndFixBfgs4, tolerance = 1e-3 )

## # with unused Hessians
## mlghFixBfgs <- maxLik( llf, gf, hf, start = startValFix, fixed = isFixed,
##    method = "BFGS" )
## expect_equal( mlgFixBfgs, mlghFixBfgs, tolerance = 1e-3 )
## mlghFixBfgs3 <- maxLik( llf, gf, hf, start = startValFix, fixed = "mu",
##    method = "BFGS" )
## expect_equal( mlghFixBfgs, mlghFixBfgs3, tolerance = 1e-3 )
## mlghFixBfgs4 <- maxLik( llf, gf, hf, start = startValFix, fixed = 1,
##    method = "BFGS" )
## expect_equal( mlghFixBfgs, mlghFixBfgs4, tolerance = 1e-3 )

## ## NM method with fixed parameters
## mlFixNm <- maxLik( llf, start = startValFix, fixed = isFixed,
##    method = "NM" )
## mlFixNm3 <- maxLik( llf, start = startValFix, fixed = "mu",
##    method = "NM" )
## expect_equal( mlFixNm, mlFixNm3, tolerance = 1e-3 )
## mlFixNm4 <- maxLik( llf, start = startValFix, fixed = 1,
##    method = "NM" )
## expect_equal( mlFixNm, mlFixNm4, tolerance = 1e-3 )
## print( mlFixNm )
## print( summary( mlFixNm ), digits = 2 )
## activePar( mlFixNm )
## AIC( mlFixNm )
## coef( mlFixNm )
## condiNumber( mlFixNm, digits = 3 )
## round( hessian( mlFixNm ), 1 )
## logLik( mlFixNm )
## maximType( mlFixNm )
## nIter( mlFixNm )
## nParam( mlFixNm )
## returnCode( mlFixNm )
## returnMessage( mlFixNm )
## round( vcov( mlFixNm ), 3 )
## logLik( summary( mlFixNm ) )
## expect_equal( mlFixBfgs[ -c(4,9,10) ], mlFixNm[ -c(4,9,10) ], tolerance = 1e-3 )
## mlIndFixNm <- maxLik( llfInd, start = startValFix, fixed = isFixed,
##    method = "NM" )
## expect_equal( mlFixNm[-4], mlIndFixNm[-c(4,12)], tolerance = 1e-3 )
## round( mlIndFixNm[[ 12 ]], 3 )
## mlIndFixNm3 <- maxLik( llfInd, start = startValFix, fixed = "mu",
##    method = "NM" )
## expect_equal( mlIndFixNm, mlIndFixNm3, tolerance = 1e-3 )
## mlIndFixNm4 <- maxLik( llfInd, start = startValFix, fixed = 1,
##    method = "NM" )
## expect_equal( mlIndFixNm, mlIndFixNm4, tolerance = 1e-3 )
## nObs( mlIndFixNm )

## # with analytical gradients
## mlgFixNm <- maxLik( llf, gf, start = startValFix, fixed = isFixed,
##    method = "NM" )
## mlgFixNm3 <- maxLik( llf, gf, start = startValFix, fixed = "mu",
##    method = "NM" )
## expect_equal( mlgFixNm, mlgFixNm3, tolerance = 1e-3 )
## mlgFixNm4 <- maxLik( llf, gf, start = startValFix, fixed = 1,
##    method = "NM" )
## expect_equal( mlgFixNm, mlgFixNm4, tolerance = 1e-3 )
## print( summary( mlgFixNm ), digits = 2 )
## expect_equal( mlFixNm, mlgFixNm, tolerance = 1e-3 )
## mlgIndFixNm <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
##    method = "NM")
## expect_equal( mlgFixNm[ ], mlgIndFixNm[ -12 ], tolerance = 1e-3 )
## round( mlgIndFixNm[[ 12 ]], 3 )

## # with unused Hessians
## mlghFixNm <- maxLik( llf, gf, hf, start = startValFix, fixed = isFixed,
##    method = "NM" )
## expect_equal( mlgFixNm, mlghFixNm, tolerance = 1e-3 )
## mlghFixNm3 <- maxLik( llf, gf, hf, start = startValFix, fixed = "mu",
##    method = "NM" )
## expect_equal( mlghFixNm, mlghFixNm3, tolerance = 1e-3 )
## mlghFixNm4 <- maxLik( llf, gf, hf, start = startValFix, fixed = 1,
##    method = "NM" )
## expect_equal( mlghFixNm, mlghFixNm4, tolerance = 1e-3 )

## ## SANN method with fixed parameters
## mlFixSann <- maxLik( llf, start = startValFix, fixed = isFixed,
##    method = "SANN" )
## mlFixSann3 <- maxLik( llf, start = startValFix, fixed = "mu",
##    method = "SANN" )
## expect_equal( mlFixSann, mlFixSann3, tolerance = 1e-3 )
## mlFixSann4 <- maxLik( llf, start = startValFix, fixed = 1,
##    method = "SANN" )
## expect_equal( mlFixSann, mlFixSann4, tolerance = 1e-3 )
## print( mlFixSann )
## print( summary( mlFixSann ), digits = 2 )
## activePar( mlFixSann )
## AIC( mlFixSann )
## coef( mlFixSann )
## condiNumber( mlFixSann, digits = 3 )
## round( hessian( mlFixSann ), 1 )
## logLik( mlFixSann )
## maximType( mlFixSann )
## nIter( mlFixSann )
## nParam( mlFixSann )
## returnCode( mlFixSann )
## returnMessage( mlFixSann )
## round( vcov( mlFixSann ), 3 )
## logLik( summary( mlFixSann ) )
## expect_equal( mlFixBfgs[ -c(4,9,10) ], mlFixSann[ -c(4,9,10) ], 
##    tolerance = 1e-3 )
## mlIndFixSann <- maxLik( llfInd, start = startValFix, fixed = isFixed,
##    method = "SANN" )
## expect_equal( mlFixSann[ ], mlIndFixSann[ -12 ], tolerance = 1e-2 )
## round( mlIndFixSann[[ 12 ]], 3 )
## nObs( mlIndFixSann )

## # with analytical gradients
## mlgFixSann <- maxLik( llf, gf, start = startValFix, fixed = isFixed,
##    method = "SANN" )
## print( summary( mlgFixSann ), digits = 2 )
## expect_equal( mlFixSann[-4], mlgFixSann[-4], tolerance = 1e-3 )
## mlgIndFixSann <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
##    method = "SANN")
## expect_equal( mlgFixSann[ ], mlgIndFixSann[ -12 ], tolerance = 1e-3 )
## round( mlgIndFixSann[[ 12 ]], 3 )

## # with unused Hessians
## mlghFixSann <- maxLik( llf, gf, hf, start = startValFix, fixed = isFixed,
##    method = "SANN" )
## expect_equal( mlgFixSann, mlghFixSann, tolerance = 1e-3 )


## ############### inequality constraints ###############
## A <- matrix( -1, nrow = 1, ncol = 2 )
## inEq <- list( ineqA = A, ineqB = 2.5 )

## ## BFGS method with inequality constraints
## mlBfgsInEq <- maxLik( llf, start = startVal, constraints = inEq,
##    method = "BFGS" )
## print( mlBfgsInEq )
## print( summary( mlBfgsInEq ), digits = 2 )
## activePar( mlBfgsInEq )
## AIC( mlBfgsInEq )
## coef( mlBfgsInEq )
## condiNumber( mlBfgsInEq, digits = 3 )
## round( hessian( mlBfgsInEq ), 1 )
## logLik( mlBfgsInEq )
## maximType( mlBfgsInEq )
## nIter( mlBfgsInEq )
## nParam( mlBfgsInEq )
## returnCode( mlBfgsInEq )
## returnMessage( mlBfgsInEq )
## round( vcov( mlBfgsInEq ), 3 )
## logLik( summary( mlBfgsInEq ) )
## mlBfgsInEqInd <- maxLik( llfInd, start = startVal, constraints = inEq,
##    method = "BFGS" )
## print( summary( mlBfgsInEqInd ), digits = 2 )
## expect_equal( mlBfgsInEq[ ], mlBfgsInEqInd[ -12 ], tolerance = 1e-3 )
## round( mlBfgsInEqInd[[ 12 ]], 3 )
## nObs( mlBfgsInEqInd )

## # with analytical gradients
## mlgBfgsInEq <- maxLik( llf, gf, start = startVal, constraints = inEq,
##    method = "BFGS" )
## expect_equal( mlBfgsInEq, mlgBfgsInEq, tolerance = 1e-3 )
## mlgBfgsInEqInd <- maxLik( llfInd, gfInd, start = startVal, constraints = inEq,
##    method = "BFGS" )
## expect_equal( mlgBfgsInEqInd[ -12 ], mlgBfgsInEq[ ], tolerance = 1e-3 )
## round( mlgBfgsInEqInd[[ 12 ]], 3 )
## mlgBfgsInEqInd2 <- maxLik( llf, gfInd, start = startVal, constraints = inEq,
##    method = "BFGS" )
## expect_equal( mlgBfgsInEqInd, mlgBfgsInEqInd2, tolerance = 1e-3 )

## # with unused Hessian
## mlghBfgsInEq <- maxLik( llf, gf, hf, start = startVal, constraints = inEq,
##    method = "BFGS" )
## expect_equal( mlgBfgsInEq, mlghBfgsInEq, tolerance = 1e-3 )

## ## NM method with inequality constraints
## mlNmInEq <- maxLik( llf, start = startVal, constraints = inEq, method = "NM" )
## print( mlNmInEq )
## print( summary( mlNmInEq ), digits = 2 )
## activePar( mlNmInEq )
## AIC( mlNmInEq )
## coef( mlNmInEq )
## condiNumber( mlNmInEq, digits = 3 )
## round( hessian( mlNmInEq ), 1 )
## logLik( mlNmInEq )
## maximType( mlNmInEq )
## nIter( mlNmInEq )
## nParam( mlNmInEq )
## returnCode( mlNmInEq )
## returnMessage( mlNmInEq )
## round( vcov( mlNmInEq ), 3 )
## logLik( summary( mlNmInEq ) )
## expect_equal( mlBfgsInEq[-c(9,10,11)], mlNmInEq[-c(9,10,11)], tolerance = 1e-3 )
## mlNmInEqInd <- maxLik( llfInd, start = startVal, constraints = inEq,
##    method = "NM" )
## print( summary( mlNmInEqInd ), digits = 2 )
## expect_equal( mlNmInEq[-4], mlNmInEqInd[-c(4,12)], tolerance = 1e-3 )
## round( mlNmInEqInd[[ 12 ]], 3 )
## nObs( mlNmInEqInd )

## # with unused analytical gradients
## mlgNmInEq <- maxLik( llf, gf, start = startVal, constraints = inEq,
##    method = "NM" )
## expect_equal( mlNmInEq, mlgNmInEq, tolerance = 1e-3 )

## # with unused analytical gradients and Hessians
## mlghNmInEq <- maxLik( llf, gf, hf, start = startVal, constraints = inEq,
##    method = "NM" )
## expect_equal( mlgNmInEq, mlghNmInEq, tolerance = 1e-3 )

## ## SANN method with inequality constraints
## mlSannInEq <- maxLik( llf, start = startVal, constraints = inEq,
##    method = "SANN" )
## print( mlSannInEq )
## print( summary( mlSannInEq ), digits = 2 )
## activePar( mlSannInEq )
## AIC( mlSannInEq )
## coef( mlSannInEq )
## condiNumber( mlSannInEq, digits = 3 )
## round( hessian( mlSannInEq ), 1 )
## logLik( mlSannInEq )
## maximType( mlSannInEq )
## nIter( mlSannInEq )
## nParam( mlSannInEq )
## returnCode( mlSannInEq )
## returnMessage( mlSannInEq )
## round( vcov( mlSannInEq ), 3 )
## logLik( summary( mlSannInEq ) )
## expect_equal( mlBfgsInEq[-c(2,3,4,9,10,11)], mlSannInEq[-c(2,3,4,9,10,11)], 
##    tolerance = 1e-3 )
## expect_equal( mlBfgsInEq[-c(3,4,9,10,11)], mlSannInEq[-c(3,4,9,10,11)], 
##    tolerance = 1e-2 )
## # with unused analytical gradients
## mlgSannInEq <- maxLik( llf, gf, start = startVal, constraints = inEq,
##    method = "SANN" )
## expect_equal( mlSannInEq, mlgSannInEq, tolerance = 1e-3 )

## # with a user-specified function to generate a new candidate point
## mlSannInEqCand <- maxLik( llf, start = startVal, constraints = inEq,
##    method = "SANN", cand = function(x)c(rnorm(1,x[1]),rnorm(1,x[2])) )
## print( summary( mlSannInEqCand ), digits = 2 )
## expect_equal( mlSannInEqCand[-c(2,3,4)], mlSannInEq[-c(2,3,4)], tolerance = 1e-3 )
## expect_equal( mlSannInEqCand, mlSannInEq, tolerance = 1e-1 )

## ############### equality constraints ###############
## eqCon <- list( eqA = A, eqB = 2.5 )

## # with analytical gradients as attribute
## mlGCon <- maxLik( llfGrad, start = startVal, constraints = eqCon )

## # with analytical gradients and Hessians
## mlghCon <- maxLik( llf, gf, hf, start = startVal, constraints = eqCon )
## expect_equal( mlGCon, mlghCon, tolerance = 1e-3 )

## # with analytical gradients and Hessians as attributes
## mlGHCon <- maxLik( llfGradHess, start = startVal, constraints = eqCon )
## expect_equal( mlGHCon, mlghCon, tolerance = 1e-3 )
## expect_equal( mlGHCon[-c(2,3,4,5,6,7,9,11)], mlGCon[-c(2,3,4,5,6,7,9,11)], 
##    tolerance = 1e-3 )
## expect_equal( mlGHCon[-c(5,6,7,9,11)], mlGCon[-c(5,6,7,9,11)], 
##    tolerance = 1e-1 )


## ## BHHH method with equality constraints
## mlBhhhCon <- maxLik( llfInd, start = startVal, constraints = eqCon,
##    method = "BHHH" )
## print( mlBhhhCon )
## print( summary( mlBhhhCon ), digits = 2 )
## activePar( mlBhhhCon )
## AIC( mlBhhhCon )
## coef( mlBhhhCon )
## condiNumber( mlBhhhCon, digits = 3 )
## round( hessian( mlBhhhCon ), 1 )
## logLik( mlBhhhCon )
## maximType( mlBhhhCon )
## nIter( mlBhhhCon )
## nParam( mlBhhhCon )
## returnCode( mlBhhhCon )
## returnMessage( mlBhhhCon )
## round( vcov( mlBhhhCon ), 3 )
## logLik( summary( mlBhhhCon ) )
## expect_equal( mlGCon[ -c( 5, 6, 7, 9, 10 ) ], mlBhhhCon[ -c( 5, 6, 7, 9, 10, 11 ) ],
##    tolerance = 5e-3 )
## mlBhhhCon[11]
## nObs( mlBhhhCon )

## # with analytical gradients
## mlgBhhhCon <- maxLik( llf, gfInd, start = startVal, constraints = eqCon,
##    method = "BHHH" )
## print( summary( mlgBhhhCon ), digits = 2 )
## expect_equal( mlBhhhCon[-c(2,3,4,5,6,7,9,11,12)], mlgBhhhCon[-c(2,3,4,5,6,7,9,11,12)],
##    tolerance = 1e-3 )
## expect_equal( mlBhhhCon[-c(5,6,7,9,12)], mlgBhhhCon[-c(5,6,7,9,12)],
##    tolerance = 1e-1 )
## mlgBhhhConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon,
##    method = "BHHH" )
## expect_equal( mlgBhhhCon, mlgBhhhConInd, tolerance = 1e-3 )

## # with analytical gradients as attribute
## mlGBhhhCon <- maxLik( llfGradInd, start = startVal, constraints = eqCon,
##    method = "BHHH" )
## print( summary( mlGBhhhCon ), digits = 2 )
## expect_equal( mlGBhhhCon, mlgBhhhCon, tolerance = 1e-3 )
## expect_equal( mlGBhhhCon[-c(2,3,4,5,6,7,9,11,12)], mlBhhhCon[-c(2,3,4,5,6,7,9,11,12)],
##    tolerance = 1e-3 )
## expect_equal( mlGBhhhCon[-c(5,6,7,9,12)], mlBhhhCon[-c(5,6,7,9,12)],
##    tolerance = 1e-1 )

## # with analytical gradients and unused Hessians
## mlghBhhhCon <- maxLik( llf, gfInd, hf, start = startVal, constraints = eqCon,
##    method = "BHHH" )
## expect_equal( mlgBhhhCon, mlghBhhhCon, tolerance = 1e-3 )

## # with analytical gradients and unused Hessians as attributes
## mlGHBhhhCon <- maxLik( llfGradHessInd, start = startVal, constraints = eqCon,
##    method = "BHHH" )
## expect_equal( mlGHBhhhCon, mlghBhhhCon, tolerance = 1e-3 )
## expect_equal( mlGHBhhhCon, mlGBhhhCon, tolerance = 1e-3 )


## ## BFGS method with equality constraints
## mlBfgsCon <- maxLik( llf, start = startVal, constraints = eqCon,
##    method = "BFGS" )
## print( mlBfgsCon )
## print( summary( mlBfgsCon ), digits = 2 )
## activePar( mlBfgsCon )
## AIC( mlBfgsCon )
## coef( mlBfgsCon )
## condiNumber( mlBfgsCon, digits = 3 )
## round( hessian( mlBfgsCon ), 1 )
## logLik( mlBfgsCon )
## maximType( mlBfgsCon )
## nIter( mlBfgsCon )
## nParam( mlBfgsCon )
## returnCode( mlBfgsCon )
## returnMessage( mlBfgsCon )
## round( vcov( mlBfgsCon ), 3 )
## logLik( summary( mlBfgsCon ) )
## mlBfgsConInd <- maxLik( llfInd, start = startVal, constraints = eqCon,
##    method = "BFGS" )
## print( summary( mlBfgsConInd ), digits = 2 )
## expect_equal( mlBfgsCon[-c(4,9)], mlBfgsConInd[-c(4,9,12)], tolerance = 1e-3 )
## mlBfgsConInd[12]
## nObs( mlBfgsConInd )

## # with analytical gradients
## mlgBfgsCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
##    method = "BFGS" )
## print( summary( mlgBfgsCon ), digits = 2 )
## expect_equal( mlBfgsCon[-c(3,4,9,11)], mlgBfgsCon[-c(3,4,9,11)], tolerance = 1e-2 )
## mlgBfgsConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon,
##    method = "BFGS" )
## expect_equal( mlgBfgsCon[], mlgBfgsConInd[-12], tolerance = 1e-3 )
## mlgBfgsConInd[12]

## # with analytical gradients and unused Hessians
## mlghBfgsCon <- maxLik( llf, gf, hf, start = startVal, constraints = eqCon,
##    method = "BFGS" )
## expect_equal( mlgBfgsCon, mlghBfgsCon, tolerance = 1e-3 )

## ## NM method with equality constraints
## mlNmCon <- maxLik( llf, start = startVal, constraints = eqCon, method = "NM", SUMTTol=0)
## print( mlNmCon )
## print( summary( mlNmCon ), digits = 2 )
## activePar( mlNmCon )
## AIC( mlNmCon )
## coef( mlNmCon )
## condiNumber( mlNmCon, digits = 3 )
## round( hessian( mlNmCon ), 1 )
## logLik( mlNmCon )
## maximType( mlNmCon )
## nIter( mlNmCon )
## nParam( mlNmCon )
## returnCode( mlNmCon )
## returnMessage( mlNmCon )
## round( vcov( mlNmCon ), 3 )
## logLik( summary( mlNmCon ) )
## mlNmConInd <- maxLik( llfInd, start = startVal, constraints = eqCon,
##    method = "NM", SUMTTol=0)
## print( summary( mlNmConInd ), digits = 2 )
## expect_equal( mlNmCon[], mlNmConInd[-12], tolerance = 1e-3 )
## mlNmConInd[12]
## nObs( mlNmConInd )

## # with unused analytical gradients
## mlgNmCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
##    method = "NM", SUMTTol=0)
## expect_equal( mlNmCon, mlgNmCon, tolerance = 1e-3 )
## mlgNmConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon,
##    method = "NM", SUMTTol=0)
## expect_equal( mlgNmCon[], mlgNmConInd[-12], tolerance = 1e-3 )
## mlgNmConInd[12]

## # with unused analytical gradients and Hessians
## mlghNmCon <- maxLik( llf, gf, hf, start = startVal, constraints = eqCon,
##    method = "NM", SUMTTol=0)
## expect_equal( mlgNmCon, mlghNmCon, tolerance = 1e-3 )

## ## SANN method with equality constraints
## mlSannCon <- maxLik( llf, start = startVal, constraints = eqCon,
##    method = "SANN", SUMTTol=0)
## print( mlSannCon )
## print( summary( mlSannCon ), digits = 2 )
## activePar( mlSannCon )
## AIC( mlSannCon )
## coef( mlSannCon )
## condiNumber( mlSannCon, digits = 3 )
## round( hessian( mlSannCon ), 1 )
## logLik( mlSannCon )
## maximType( mlSannCon )
## nIter( mlSannCon )
## nParam( mlSannCon )
## returnCode( mlSannCon )
## returnMessage( mlSannCon )
## round( vcov( mlSannCon ), 3 )
## logLik( summary( mlSannCon ) )
## expect_equal( mlSannCon[ -c(2,3,4,5,6,9,10,11) ], mlBfgsCon[ -c(2,3,4,5,6,9,10,11) ],
##    tolerance = 1e-3 )
## expect_equal( mlSannCon[ -c(3,4,5,6,9,10,11) ], mlBfgsCon[ -c(3,4,5,6,9,10,11) ],
##    tolerance = 1e-2 )

## # with unused analytical gradients
## mlgSannCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
##    method = "SANN", SUMTTol=0)
## expect_equal( mlSannCon, mlgSannCon, tolerance = 1e-3 )

## # with a user-specified function to generate a new candidate point
## mlSannConCand <- maxLik( llf, start = startVal, constraints = eqCon,
##    method = "SANN", cand = function(x)c(rnorm(1,x[1]),rnorm(1,x[2])) )
## print( summary( mlSannConCand ), digits = 2 )
## expect_equal( mlSannConCand[-c(1,2,3,4,11)], mlSannCon[-c(1,2,3,4,11)], 
##    tolerance = 1e-3 )
## expect_equal( mlSannConCand[-c(2,3,4,11)], mlSannCon[-c(2,3,4,11)], 
##    tolerance = 1e-1 )


## ## test for method "estfun"
## library( sandwich )
## try( estfun( ml ) )
## estfun( mlInd )[ 1:5, ]
## estfun( mlgInd )[ 1:5, ]
## estfun( mlBHHH )[ 1:5, ]
## estfun( mlgBHHH )[ 1:5, ]
## estfun( mlIndBFGS )[ 1:5, ]
## estfun( mlgIndBFGS )[ 1:5, ]
## estfun( mlIndNM )[ 1:5, ]
## estfun( mlgIndNM )[ 1:5, ]
## estfun( mlIndSANN )[ 1:5, ]
## estfun( mlgIndSANN )[ 1:5, ]
## estfun( mlIndFix )[ 1:5, ]
## estfun( mlgIndFix )[ 1:5, ]
## estfun( mlFixBHHH )[ 1:5, ]
## estfun( mlgFixBHHH )[ 1:5, ]
## estfun( mlIndFixBfgs )[ 1:5, ]
## estfun( mlgIndFixBfgs )[ 1:5, ]
## estfun( mlIndFixNm )[ 1:5, ]
## estfun( mlgIndFixNm )[ 1:5, ]
## estfun( mlIndFixSann )[ 1:5, ]
## estfun( mlgIndFixSann )[ 1:5, ]
## estfun( mlBfgsInEqInd )[ 1:5, ]
## estfun( mlgBfgsInEqInd )[ 1:5, ]
## estfun( mlNmInEqInd )[ 1:5, ]
## estfun( mlBhhhCon )[ 1:5, ]
## estfun( mlgBhhhCon )[ 1:5, ]
## estfun( mlBfgsConInd )[ 1:5, ]
## estfun( mlgBfgsConInd )[ 1:5, ]
## estfun( mlNmConInd )[ 1:5, ]
## estfun( mlgNmConInd )[ 1:5, ]


## ## test for method "bread"
## try( bread( ml ) )
## round( bread( mlInd ), 2 )
## round( bread( mlgInd ), 2 )
## round( bread( mlBHHH ), 2 )
## round( bread( mlgBHHH ), 2 )
## round( bread( mlIndBFGS ), 2 )
## round( bread( mlgIndBFGS ), 2 )
## round( bread( mlIndNM ), 2 )
## round( bread( mlgIndNM ), 2 )
## round( bread( mlIndSANN ), 2 )
## round( bread( mlgIndSANN ), 2 )
## round( bread( mlIndFix ), 2 )
## round( bread( mlgIndFix ), 2 )
## round( bread( mlFixBHHH ), 2 )
## round( bread( mlgFixBHHH ), 2 )
## round( bread( mlIndFixBfgs ), 2 )
## round( bread( mlgIndFixBfgs ), 2 )
## round( bread( mlIndFixNm ), 2 )
## round( bread( mlgIndFixNm ), 2 )
## round( bread( mlIndFixSann ), 2 )
## round( bread( mlgIndFixSann ), 2 )
## round( bread( mlBfgsInEqInd ), 2 )
## round( bread( mlgBfgsInEqInd ), 2 )
## round( bread( mlNmInEqInd ), 2 )
## round( bread( mlBhhhCon ), 2 )
## round( bread( mlgBhhhCon ), 2 )
## round( bread( mlBfgsConInd ), 2 )
## round( bread( mlgBfgsConInd ), 2 )
## round( bread( mlNmConInd ), 2 )
## round( bread( mlgNmConInd ), 2 )


## ## test for method "sandwich"
## try( sandwich( ml ) )
## printSandwich <- function( x ) {
##    print( round( sandwich( x ), 2 ) )
##    tmp <- expect_equal( sandwich( x ), vcov( x ) )
##    if( isTRUE( tmp ) ) {
##       print( tmp )
##    }
## }
## printSandwich( mlInd )
## printSandwich( mlgInd )
## printSandwich( mlBHHH )
## printSandwich( mlgBHHH )
## printSandwich( mlIndBFGS )
## printSandwich( mlgIndBFGS )
## printSandwich( mlIndNM )
## printSandwich( mlgIndNM )
## printSandwich( mlIndSANN )
## printSandwich( mlgIndSANN )
## printSandwich( mlIndFix )
## printSandwich( mlgIndFix )
## printSandwich( mlFixBHHH )
## printSandwich( mlgFixBHHH )
## printSandwich( mlIndFixBfgs )
## printSandwich( mlgIndFixBfgs )
## printSandwich( mlIndFixNm )
## printSandwich( mlgIndFixNm )
## printSandwich( mlIndFixSann )
## printSandwich( mlgIndFixSann )
## printSandwich( mlBfgsInEqInd )
## printSandwich( mlgBfgsInEqInd )
## printSandwich( mlNmInEqInd )
## printSandwich( mlBhhhCon )
## printSandwich( mlgBhhhCon )
## printSandwich( mlBfgsConInd )
## printSandwich( mlgBfgsConInd )
## printSandwich( mlNmConInd )
## printSandwich( mlgNmConInd )

## ### ---------- convergence tolerance parameters ----------
## a <- maxNR(llf, gf, hf, start=startVal, tol=1e-3, reltol=0, gradtol=0, iterlim=10)
## returnMessage(a)  # should stop with code 2: tolerance
## a <- maxNR(llf, gf, hf, start=startVal, tol=0, reltol=1e-3, gradtol=0, iterlim=10)
## returnMessage(a)  # 8: relative tolerance
## a <- maxNR(llf, gf, hf, start=startVal, tol=0, reltol=0, gradtol=1e-3, iterlim=10)
## returnMessage(a)  # 1: gradient
## a <- maxNR(llf, gf, hf, start=startVal, tol=0, reltol=0, gradtol=0, iterlim=10)
## returnMessage(a)  # 4: iteration limit
