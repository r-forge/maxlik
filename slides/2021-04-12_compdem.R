library(maxLik)

### ---------- Simple example to show the basics ----------
### 
x <- rnorm(100)  # data.  true mu = 0, sigma = 1
loglik <- function(theta) {
   mu <- theta[1]
   sigma <- theta[2]
   sum(dnorm(x, mean=mu, sd=sigma, log=TRUE))
}
m <- maxLik(loglik, start=c(mu=1, sigma=2))
                           # give start value somewhat off
summary(m)
## Auxiliary functions (there are more)
coef(m)
stdEr(m)

### ---------- Demonstrate how numeric gradient fails while analytic works ----------
### Do a simple OLS with very unequal variable scale
### Use matrix approach (makes life _much_ easier)
### 
## create 3 variables with very different scale
X <- cbind(rnorm(100), rnorm(100, sd=1e3), rnorm(100, sd=1e7))
## note: correct coefficients are 1, 1, 1
y <- X %*% c(1,1,1) + rnorm(100)
negSSE <- function(beta) {
   e <- y - X %*% beta
   -crossprod(e)
                           # note '-': we are maximizing
}
m <- maxLik(negSSE, start=c(0,0,0))
                           # give start values a bit off
summary(m, eigentol=1e-15)  # total mess

## Now code analytic gradient
grad <- function(beta) {
   2*t(y - X %*% beta) %*% X
}
m <- maxLik(negSSE, grad=grad, start=c(0,0,0))
summary(m, eigentol=1e-15)  # works :-)

## Add analytic hessian
hess <- function(beta) {
   -2*crossprod(X)
}
m <- maxLik(negSSE, grad=grad, hess=hess, start=c(0,0,0))
summary(m, eigentol=1e-15)  # works too

## Demonstrate gradient, Hessian as attributes
negSSEA <- function(beta) {
   ## negative SSE with attributes
   e <- y - X %*% beta  # we will re-use 'e'
   sse <- -crossprod(e)
                           # note '-': we are maximizing
   attr(sse, "gradient") <- 2*t(e) %*% X
   attr(sse, "Hessian") <- -2*crossprod(X)
   sse
}
m <- maxLik(negSSEA, start=c(0,0,0))
summary(m, eigentol=1e-15)

### ---------- comparing analytic, numeric derivatives ----------
### 
compareDerivatives(negSSE, grad, t0=c(0,0,0))
                           # 't0' is the parameter value

### ---------- other optimizers ----------
### 
m <- maxLik(loglik, start=c(mu=1, sigma=2),
            method="BFGS")
summary(m)

### ---------- BHHH ----------
###
## log-likelihood of normal distribution
loglik <- function(theta) {
   mu <- theta[1]
   sigma <- theta[2]
   N <- length(x)
   -N*log(sqrt(2*pi)) - N*log(sigma) - sum(0.5*(x - mu)^2/sigma^2)
                           # sum over observations
}
## gradient coded with observation-wise values
gradlikB <- function(theta) {
   ## BHHH-compatible gradient
   mu <- theta[1]
   sigma <- theta[2]
   N <- length(x)  # number of observations
   gradient <- matrix(0, N, 2)  # gradient is matrix:
                           # N datapoints (rows), 2 components
   gradient[, 1] <- (x - mu)/sigma^2
                           # first column: derivative wrt mu
   gradient[, 2] <- -1/sigma + (x - mu)^2/sigma^3
                           # second column: derivative wrt sigma
   gradient
}
m <- maxLik(loglik, gradlikB, start=c(mu=1, sigma=2),
            method="BHHH")
summary(m)

### ---------- compute numeric gradient observation-wise ----------
### 
loglikB <- function(theta) {
   mu <- theta[1]
   sigma <- theta[2]
   -log(sqrt(2*pi)) - log(sigma) - 0.5*(x - mu)^2/sigma^2
                           # no summing here
                           # also no 'N*' terms as we work by
                           # individual observations
}
m <- maxLik(loglikB, start=c(mu=1, sigma=2),
            method="BHHH")
summary(m)

### ---------- control options ----------
###
## print a lot, only 3 iterations
m <- maxLik(loglikB, start=c(mu=1, sigma=2),
            method="BHHH",
            control=list(printLevel=3, iterlim=2))
summary(m)

## stop only at absolute tolerance
m <- maxLik(loglikB, start=c(mu=1, sigma=2),
            method="BHHH",
            control=list(reltol=0, gradtol=0))
summary(m)

### ---------- additional arguments for log-likelihood ----------
### 
loglik <- function(theta, x) {
   mu <- theta[1]
   sigma <- theta[2]
   sum(dnorm(x, mean=mu, sd=sigma, log=TRUE))
}
m <- maxLik(loglik, start=c(mu=1, sigma=2), x=x)
                           # named argument 'x' will be passed
                           # to loglik
summary(m)

### ---------- optimize non-log likelihood ----------
### 
f <- function(theta) {
   x <- theta[1]
   y <- theta[2]
   exp(-x^2 - y^2)
                           # optimum at (0, 0)
}
m <- maxBFGS(f, start=c(1,1))
                           # give start value a bit off
summary(m)

### ---------- testing condition numbers ----------
### 
## create 3 variables, two independent, third collinear
x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- x1 + x2 + rnorm(100, sd=1e-6)  # highly correlated w/x1, x2
X <- cbind(x1, x2, x3)
y <- X %*% c(1, 1, 1) + rnorm(100)
m <- maxLik(negSSEA, start=c(x1=0, x2=0, x3=0))
                           # negSSEA: negative sum of squared errors
                           # with gradient, hessian attribute
summary(m)
condiNumber(X)  # third variable crap

## linearly separated logistic regression case
x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- rnorm(100)
X <- cbind(x1, x2, x3)
y <- X %*% c(1, 1, 1) > 0
                           # y values 1/0 linearly separated
loglik <- function(beta) {
   link <- X %*% beta
   sum(ifelse(y > 0, plogis(link, log=TRUE),
              plogis(-link, log=TRUE)))
}
m <- maxLik(loglik, start=c(x1=0, x2=0, x3=0))
summary(m)
condiNumber(X)  # looks good
condiNumber(hessian(m))  # third variable crap

### ---------- fixed parameters ----------
### 
x <- rnorm(100)
loglik <- function(theta) {
   mu <- theta[1]
   sigma <- theta[2]
   sum(dnorm(x, mean=mu, sd=sigma, log=TRUE))
}
m <- maxLik(loglik, start=c(mu=1, sigma=1),
            fixed="sigma")
                           # fix the component named 'sigma'
summary(m)

### ---------- contrained maximization ----------
### 
f <- function(theta) {
   x <- theta[1]
   y <- theta[2]
   exp(-x^2 - y^2)
                           # optimum at (0, 0)
}
## equality constraints
A <- matrix(c(1, 1), ncol=2)
B <- -1
m <- maxNR(f, start=c(1,1),
           constraints=list(eqA=A, eqB=B))
summary(m)

## inequality constraints
A <- matrix(c(1, 1), ncol=2)
B <- -1
m <- maxBFGS(f, start=c(1,1),
             constraints=list(ineqA=A, ineqB=B))
summary(m)

## multiple inequality constraints
A <- matrix(c(1, 1, 1, -1), ncol=2)
B <- c(-1, -1)
m <- maxBFGS(f, start=c(2, 0),
             constraints=list(ineqA=A, ineqB=B))
summary(m)
