\name{vcov.maxLik}
\alias{vcov.maxLik}

\title{Variance Covariance Matrix of maxLik objects}
\description{
   Extract variance-covariance matrices of objects of class \code{\link{maxLik}}.
}
\usage{
   \method{vcov}{maxLik}( object, ... )
}
\arguments{
   \item{object}{an object of class \code{probit} or \code{maxLik}.}
   \item{\dots}{further arguments (currently ignored).}
}
\value{
   the estimated variance covariance matrix of the coefficients.
}
\author{
   Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de},
   Ott Toomet \email{otoomet@ut.ee}
}

\seealso{\code{\link[stat]{vcov}}, \code{\link{maxLik}}.}
\keyword{methods}
