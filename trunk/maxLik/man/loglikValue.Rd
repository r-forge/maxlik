\name{loglikValue}
\alias{loglikValue}
\alias{loglikValue.default}
\alias{loglikValue.maxLik}
\title{Return the log likelihood value}
\description{
  Return the log likelihood value.  Objects should be models,
  estimated by Maksimum Likelihood, or other statistical models where
  the likelihood makes sense.
}
\usage{
loglikValue(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object, usually a model estimated with Maximum
    Likelihood}
  \item{...}{additional arguments to methods}
}
\details{
  The default method returns the component \code{loglik} 
}
\value{
  A single numeric, log likelihood of the estimated model
}
\author{Ott Toomet \email{otoomet@ut.ee}}
\seealso{\code{\link{maxLik}}}
\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
## Estimate with analytic gradient and hessian
a <- maxLik(loglik, gradlik, hesslik, start=1)
## print log likelihood value
loglikValue( a )
}
\keyword{methods}
