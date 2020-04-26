\name{storedValues}
\alias{storedValues}
\alias{storedValues.maxim}
\title{Return the stored values of optimization} 
\description{
  Retrieve the objective function value for each iteration if stored
  during the optimization.
}
\usage{
storedValues(x, \dots)
\method{storedValues}{maxim}(x, \dots)
}
\arguments{
  \item{x}{a result of maximization, created
    by \code{\link{maxLik}}, \code{\link{maxNR}} or another optimizer.}
  \item{\dots}{further arguments for other methods}
}
\details{
  This is a generic method.
  If asked by control parameter \code{storeValues=TRUE}, certain
  optimization methods store the objective function value.  This method
  retrieves the stored values.
}
\value{
  a numeric vector, one value for each iteration.
}
\author{Ott Toomet}
\seealso{\code{\link{maxNR}},
  \code{\link{maxControl}}
  }
\examples{
## Estimate the exponential distribution parameter
t <- rexp(100, 2)
loglik <- function(theta, index) sum(log(theta) - theta*t[index])
## Estimate with numeric gradient and numeric Hessian
a <- maxSGA(loglik, start=1, control=list(storeValues=TRUE), nObs=100)
storedValues(a)
}
\keyword{methods}
