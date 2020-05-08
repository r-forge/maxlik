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
    log-likelihood functions and it requires the gradient/log-likelihood to
    be supplied by individual observations.

    Hessian matrix is not often used for optimization problems where one
    applies SGA, but even if one is not interested in standard errors,
    it may provide useful information about the model performance.  If
    computed by finite-difference method, the Hessian computation may be
    very slow.
  }
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
      \item{SGA_clip}{\code{NULL}, gradient clipping threshold.  This is
	the max allowed squared Euclidean norm of the gradient.  If the
	actual norm of the gradient exceeds (square root of) this
	threshold, the gradient will be scaled back accordingly while
	preserving its direction.  \code{NULL} means no clipping.
      }
      \item{SGA_momentum}{0, numeric momentum parameter for SGA.  Must lie
	in interval \eqn{[0,1]}{[0,1]}.  See details.
      }
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

    If \code{SGA_momentum} is positive, the algorithm updates the parameters
    \eqn{\theta}{theta} in two steps.  First, the momentum is used to update
    the \dQuote{velocity} \eqn{v}{v} as
    \eqn{v \leftarrow \mathrm{momentum}\cdot v + \mathrm{learning
	rate}\cdot \nabla f(\theta)}{v <- momentum*v + learning
      rate* gradient f(theta)}, and thereafter the parameter
    \eqn{\theta}{theta} is updates as
    \eqn{\theta \leftarrow \theta + v}{theta <- theta + v}.  Initial
    velocity is set to 0.
    
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
  \item{\code{\link{coef}}}{estimated parameter value.}
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

\references{
    Bottou, L.; Curtis, F. & Nocedal, J.:
    Optimization Methods for
    Large-Scale Machine Learning \emph{SIAM Review}, 2018, \bold{60}, 223--311.

    Henningsen, A. and Toomet, O. (2011): maxLik: A package for maximum likelihood
    estimation in R \emph{Computational Statistics} \bold{26}, 443--458
}

\author{Ott Toomet, Arne Henningsen,
   function \code{maxBFGSR} was originally developed by Yves Croissant
   (and placed in 'mlogit' package)
 }

 \seealso{
   \code{\link{maxNR}} for Newton-Raphson, a popular Hessian-based maximization;
   \code{\link{maxBFGS}} for maximization using the BFGS, Nelder-Mead (NM),
   and Simulated Annealing (SANN) method (based on \code{\link{optim}}),
   also supporting inequality constraints;
   \code{\link{maxLik}} for a general framework for maximum likelihood
   estimation (MLE);
   \code{\link{optim}} for different gradient-based optimization
   methods.
 }

 \examples{
## estimate the exponential distribution parameter by ML
set.seed(1)
t <- rexp(100, 2)
loglik <- function(theta, index) sum(log(theta) - theta*t[index])
## Note the log-likelihood and gradient are summed over observations
gradlik <- function(theta, index) sum(1/theta - t[index])
## Estimate with finite-difference gradient and Hessian
a <- maxSGA(loglik, gradlik, start=1, control=list(iterlim=1000,
            SGA_batchSize=10), nObs=100)
            # note that loglik is not really needed, and is not used
            # here, unless more print verbosity is asked
summary(a)
##
## demonstrate the usage of index, and using
## fn for computing the objective function on validation data.
## create a linear model where variables are very unequally scaled
##
## OLS loglik function: compute the function value on validation data only
loglik <- function(beta, index) {
   e <- yValid - XValid \%*\% beta
   -crossprod(e)/length(y)
}
## OLS gradient: compute it on training data only
gradlik <- function(beta, index)  {
   e <- yTrain[index] - XTrain[index,,drop=FALSE] \%*\% beta
   l <- crossprod(e)
   g <- t(-2*t(XTrain[index,,drop=FALSE]) \%*\% e)
   -g/length(index)
}
N <- 1000
## two random variables: one with scale 1, the other with 100
X <- cbind(rnorm(N), rnorm(N, sd=100))
beta <- c(1, 1)  # true parameter values
y <- X \%*\% beta + rnorm(N, sd=0.2)
## training-validation split
iTrain <- sample(N, 0.8*N)
XTrain <- X[iTrain,,drop=FALSE]
XValid <- X[-iTrain,,drop=FALSE]
yTrain <- y[iTrain]
yValid <- y[-iTrain]
##
cat("Analytic training solution:\n")
print(solve(crossprod(XTrain)) \%*\% crossprod(XTrain, yTrain))
## do this without momentum: learning rate must stay small for the gradient not to explode
cat("  No momentum:\n")
a <- maxSGA(loglik, gradlik, start=c(10,10),
           control=list(printLevel=1, iterlim=1000,
                        SGA_batchSize=NULL, SGA_learningRate=0.0001, SGA_momentum=0,
                        storeValues=TRUE), nObs=length(yTrain))
print(summary(a))  # the first component is off, the second close to the true value
## do with momentum 0.99
cat("  Momentum 0.99:\n")
a <- maxSGA(loglik, gradlik, start=c(10,10),
           control=list(printLevel=1, iterlim=1000,
                        SGA_batchSize=NULL, SGA_learningRate=0.0001, SGA_momentum=0.99,
                        # no momentum
                        storeValues=TRUE), nObs=length(yTrain))
print(summary(a))  # close to true value
}
\keyword{optimize}
