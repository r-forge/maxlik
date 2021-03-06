\name{translogDeriv}
\alias{translogDeriv}
\title{Derivatives of a translog function}

\description{
   Calculate the derivatives of a translog function.
}

\usage{
   translogDeriv( xNames, data, allCoef, allCoefCov = NULL,
   quadHalf = TRUE, logValues = FALSE )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data.}
   \item{allCoef}{vector containing all coefficients.}
   \item{allCoefCov}{optional covariance matrix of the coefficients.}
   \item{quadHalf}{logical. Multiply the quadratic terms by one half?}
   \item{logValues}{logical. Are the values in \code{data} already logged?}
}

\value{
   a list of class \code{translogDeriv} containing following objects:
   \item{deriv}{data frame containing the derivatives.}
   \item{variance}{data frame containing the variances of the derivatives
      (not implemented yet).}
   \item{stdDev}{data frame containing the standard deviations of the derivatives
      (not implemented yet).}
}

\seealso{
   \code{\link{translogEst}}, \code{\link{translogCalc}} nad
   \code{\link{translogHessian}}
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of variable inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
   # a time trend to account for technical progress:
   germanFarms$time <- c(1:20)

   # estimate a quadratic production function
   estResult <- translogEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms )

   # compute the marginal products of the inputs
   margProducts <- translogDeriv( c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms, estResult$allCoef, estResult$allCoefCov )
   margProducts$deriv
}

\keyword{models}
