\name{maxSGA}
\alias{maxSGA}
\title{Stochastic Gradient Ascent}
\description{
  Stochastic Gradient Ascent and related methods
}
\usage{
maxSGA(fn, grad = NULL, hess = NULL, start,
nObs,
      constraints = NULL, finalHessian = FALSE, 
      fixed = NULL, activePar = NULL, control=NULL, ... )
}
\arguments{
  \item{fn}{the function to be maximized.
    It must have the parameter vector as the first argument, and it must
    have an argument \code{index} to specify the set of selected
    observations. 
    It must return either a single number, or a numeric vector (this is
    is summed internally).
    If the parameters are out of range, \code{fn} should
    return \code{NA}.  See details for constant parameters.

    \code{fn} may also return attributes "gradient" and/or "hessian".
    If these attributes are set, the algorithm uses the corresponding
    values as
    gradient and Hessian.
  }
  \item{grad}{gradient of the objective function.
    It must have the parameter vector as the first argument, and it must
    have an argument \code{index} to specify the set of selected
    observations.
    It must return either a gradient vector of the objective function,
    or a matrix, where \emph{columns} correspond to individual parameters.
    The column sums are treated as gradient components.
    If \code{NULL}, finite-difference gradients are computed.
    If BHHH method is used, \code{grad} must return a matrix,
    where rows corresponds to the gradient vectors for individual
    observations and the columns to the individual parameters.
    If \code{fn} returns an object with attribute \code{gradient},
    this argument is ignored.

    If \code{grad} is not supplied, it is derived by finite-difference
    method from \code{fn}.  However, this is only adiseable for
    small-scale tests, not for any production run.  Obviously, \code{fn}
    must be correctly defined in that case.
  }
  \item{hess}{Hessian matrix of the function.  Mainly for compatibility
    reasons, only used for computing the final Hessian if asked to do
    so by setting \code{finalHessian} to \code{TRUE}. 
    It must have the parameter vector as the first argument and
    it must return the Hessian matrix of the objective function.
    If missing, finite-difference Hessian, based on \code{gradient},
    is computed.
  }
  \item{start}{initial parameter values.  If start values
    are named, those names are also carried over to the results.}
  \item{nObs}{number of observations.  This is used to select the
    random individual batches of observations.  The resulting batch
    indices are forwarded to the \code{grad} function through the
    arument \code{index}.}
  \item{constraints}{either \code{NULL} for unconstrained optimization
    or a list with two components.  The components may be either
    \code{eqA} and \code{eqB} for equality-constrained optimization
    \eqn{A \theta + B = 0}{A \%*\% theta + B = 0}; or \code{ineqA} and
    \code{ineqB} for inequality constraints \eqn{A \theta + B > 0}{A
      \%*\% theta + B > 0}.  More
       than one
       row in \code{ineqA} and \code{ineqB} corresponds to more than
       one linear constraint, in that case all these must be zero
       (equality) or positive (inequality constraints).
       The equality-constrained problem is forwarded
    to \code{\link{sumt}}, the inequality-constrained case to
    \code{\link{constrOptim2}}.
  }
  \item{finalHessian}{how (and if) to calculate the final Hessian.  Either
    \code{FALSE} (do not calculate), \code{TRUE} (use analytic/finite-difference
    Hessian) or \code{"bhhh"}/\code{"BHHH"} for the information equality
    approach.  The latter approach is only suitable for maximizing
    log-likelihood functions.  It requires the gradient/log-likelihood to
    be supplied by individual observations.
    Note that computing the (actual, not BHHH) final Hessian
    does not carry any extra penalty for the NR method,
    but does for the other methods.}
  \item{fixed}{parameters to be treated as constants at their
    \code{start} values.  If present, it is treated as an index vector of
    \code{start} parameters.}
  \item{activePar}{this argument is retained for backward compatibility only;
    please use argument \code{fixed} instead.}
  \item{control}{list of control parameters.  The control parameters
    used by these optimizers are
    \describe{
      \item{gradtol}{stopping condition.  Stop if norm of the gradient is
	less than \code{gradtol}.  Default 0, i.e. do not use this
	condition.  Return code 1.  This condition is useful if the
	objective is to drive full batch gradient to zero on training data.
	However, it is not a good objective in case of the stochastic
	gradient, and if the objective is to optimize the objective on
	validation data.
      }
      \item{SGA_learningRate}{step size the SGA algorithm takes in the
	gradient direction.  If 1, it steps by the gradient value.  A
	good value is typically 0.01--0.3}
      \item{SGA_batchSize}{SGA batch size, a positive
	integer.  If \code{NULL}, full batch is taken.}
      \item{iterlim}{stopping condition.  Stop if more than \code{iterlim}
	iterations, return \code{code=4}.
	Epoch is a set of iterations that cycles through all
	observations.  In case of full batch, iterations and epochs are
	equivalent. 
      }
      \item{printLevel}{this argument determines the level of
	printing which is done during the optimization process. The default
	value 0 means that no printing occurs, 1 prints the
	initial and final details, 2 prints all the
	main tracing information for every iteration.  Higher
	values will result in even more output.
      }
    }
  }
  \item{\dots}{further arguments to \code{fn}, \code{grad} and
    \code{hess}.
    Further arguments to \code{maxBHHH} are also passed to
    \code{maxNR}.
    To maintain compatibility with the earlier versions, \dots also passes a
    number of control options (\code{tol}, \code{reltol},
    \code{gradtol}, \code{steptol},
    \code{lambdatol},  \code{qrtol}, \code{iterlim}) to the optimizers.
    }
  }
  \details{
    Gradient Ascent (GA) is a optimization method where the algorithm
    repeatedly takes small steps in the gradient's direction.
    In case of Stochastic GA (SGA), the gradient is not computed on the
    full set of observations but a random subset, \emph{batch},
    potentially a single observation only.  This converges much faster
    than when using all observation in certain circumstances (see
    \cite{Bottou et al, 2018}).

    The function \code{fn} is not directly used for optimization, only
    for printing or eventually as a stopping condition.  In this sense
    it is up to the user to decide what (if anything the function
    returns).  Normally it is useful for \code{fn} to compute the
    objective function on either full training data, or validation data,
    and just ignore the \code{index} argument.  However, one may also
    choose to select the observations determined byt the index to
    compute the objective function on the current data batch.

    Does it support contraints?
  }

\value{
  object of class "maxim".  Data can be extracted through the following
  methods: 
  \item{\code{\link{maxValue}}}{\code{fn} value at maximum (the last calculated value
    if not converged.)}
  \item{coef}{estimated parameter value.}
  \item{gradient}{vector, last calculated gradient value.  Should be
    close to 0 in case of normal convergence.}
  \item{estfun}{matrix of gradients at parameter value \code{estimate}
      evaluated at each observation (only if \code{grad} returns a matrix
      or \code{grad} is not specified and \code{fn} returns a vector).}
  \item{hessian}{Hessian at the maximum (the last calculated value if
    not converged).}
  \item{returnCode}{return code:
    \itemize{
    \item{1}{ gradient close to zero (normal convergence).}
    \item{2}{ successive function values within tolerance limit (normal
      convergence).}
    \item{3}{ last step could not find higher value (probably not
      converged).  This is related to line search step getting too
      small, usually because hitting the boundary of the parameter
      space.  It may also be related to attempts to move to a wrong
      direction because of numerical errors.  In some cases it can be
      helped by changing \code{steptol}.}
    \item{4}{ iteration limit exceeded.}
    \item{5}{ Infinite value.}
    \item{6}{ Infinite gradient.}
    \item{7}{ Infinite Hessian.}
    \item{8}{ Successive function values withing relative tolerance
      limit (normal convergence).}
    \item{9}{ (BFGS) Hessian approximation cannot be improved because of
      gradient did not change.  May be related to numerical
      approximation problems or wrong analytic gradient.}
    \item{100}{ Initial value out of range.}
    }
  }
  \item{\code{\link{storedValues}}}{return values stored at each epoch}
  \item{returnMessage}{ a short message, describing the return code.}
  \item{activePar}{logical vector, which parameters are optimized over.
    Contains only \code{TRUE}-s if no parameters are fixed.}
  \item{nIter}{number of iterations.}
  \item{maximType}{character string, type of maximization.}
  \item{maxControl}{the optimization control parameters in the form of a
    \code{\link[maxLik:MaxControl-class]{MaxControl}} object.}

  The following components can only be extracted directly (with \code{\$}):
  \item{last.step}{a list describing the last unsuccessful step if
    \code{code=3} with following components:
    \itemize{
    \item{theta0}{ previous parameter value}
    \item{f0}{ \code{fn} value at \code{theta0}}
    \item{climb}{ the movement vector to the maximum of the quadratic approximation}
    }
  }
  \item{constraints}{A list, describing the constrained optimization
    (\code{NULL} if unconstrained).  Includes the following components:
    \itemize{
      \item{type}{ type of constrained optimization}
      \item{outer.iterations}{ number of iterations in the constraints step}
      \item{barrier.value}{ value of the barrier function}
    }
  }
}

\section{Warning}{
    No attempt is made to ensure that user-provided analytic
    gradient/Hessian is correct.  The users are
    encouraged to use \code{\link{compareDerivatives}} function,
    designed for this purpose.  If analytic gradient/Hessian are wrong,
    the algorithm may not converge, or may converge to a wrong point.

    As the BHHH method
    uses the likelihood-specific information equality,
    it is only suitable for maximizing log-likelihood functions!

    Quasi-Newton methods, including those mentioned above, do not work
    well in non-concave regions.  This is especially the case with the
    implementation in \code{maxBFGSR}.  The user is advised to
    experiment with various tolerance options to achieve convergence.
  }
  \references{
    Bottou, L.; Curtis, F. & Nocedal, J.:
    Optimization Methods for
    Large-Scale Machine Learning \emph{SIAM Review}, 2018, \bold{60}, 223--311.

    Henningsen, A. and Toomet, O. (2011): maxLik: A package for maximum likelihood
    estimation in R \emph{Computational Statistics} \bold{26}, 443--458
}

\author{Ott Toomet, Arne Henningsen,
   function \code{maxBFGSR} was originally developed by Yves Croissant
   (and placed in 'mlogit' package)}

\seealso{\code{\link{maxLik}} for a general framework for maximum likelihood
   estimation (MLE);
   \code{\link{maxBHHH}} for maximizations using the Berndt, Hall, Hall,
   Hausman (1974) algorithm (which is a wrapper function to \code{maxNR});
   \code{\link{maxBFGS}} for maximization using the BFGS, Nelder-Mead (NM),
   and Simulated Annealing (SANN) method (based on \code{\link{optim}}),
   also supporting inequality constraints;
   \code{\link{nlm}} for Newton-Raphson optimization; and
  \code{\link{optim}} for different gradient-based optimization
  methods.}

\examples{
## estimate the exponential distribution parameter by ML
t <- rexp(100, 2)
loglik <- function(theta) sum(log(theta) - theta*t)
## Note the log-likelihood and gradient are summed over observations
gradlik <- function(theta) sum(1/theta - t)
hesslik <- function(theta) -100/theta^2
## Estimate with finite-difference gradient and Hessian
a <- maxNR(loglik, start=1, control=list(printLevel=2))
summary(a)
## You would probably prefer 1/mean(t) instead ;-)
## Estimate with analytic gradient and Hessian
a <- maxNR(loglik, gradlik, hesslik, start=1)
summary(a)

## BFGS estimation with finite-difference gradient
a <- maxBFGSR( loglik, start=1 )
summary(a)

## For the BHHH method we need likelihood values and gradients
## of individual observations
loglikInd <- function(theta) log(theta) - theta*t
gradlikInd <- function(theta) 1/theta - t
## Estimate with analytic gradient
a <- maxBHHH(loglikInd, gradlikInd, start=1)
summary(a)

##
## Example with a vector argument:  Estimate the mean and
## variance of a random normal sample by maximum likelihood
## Note: you might want to use maxLik instead
##
loglik <- function(param) {
  mu <- param[1]
  sigma <- param[2]
  ll <- -0.5*N*log(2*pi) - N*log(sigma) - sum(0.5*(x - mu)^2/sigma^2)
  ll
}
x <- rnorm(100, 1, 2) # use mean=1, stdd=2
N <- length(x)
res <- maxNR(loglik, start=c(0,1)) # use 'wrong' start values
summary(res)
##
## The previous example with named parameters and fixed values
##
resFix <- maxNR(loglik, start=c(mu=0, sigma=1), fixed="sigma")
summary(resFix)  # 'sigma' is exactly 1.000 now.
###
### Constrained optimization
###
## We maximize exp(-x^2 - y^2) where x+y = 1
hatf <- function(theta) {
  x <- theta[1]
  y <- theta[2]
  exp(-(x^2 + y^2))
  ## Note: you may prefer exp(- theta \%*\% theta) instead
}
## use constraints: x + y = 1
A <- matrix(c(1, 1), 1, 2)
B <- -1
res <- maxNR(hatf, start=c(0,0), constraints=list(eqA=A, eqB=B),
             control=list(printLevel=1))
print(summary(res))
}
\keyword{optimize}
