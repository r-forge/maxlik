### tests for stochastic gradient ascent
### Test the following things:
###
### 1. basic SGA
### 2. SGA full batch
### 3. SGA, no gradient supplied

library(maxLik)
library(testthat)

## ---------- OLS 
## log-likelihood function(s):
## return log-likelihood on validation data
loglik <- function(beta, index) {
   e <- yValid - XValid %*% beta
   -crossprod(e)/length(y)
}
## gradlik: work on training data
gradlik <- function(beta, index) {
   e <- yTrain[index] - XTrain[index,,drop=FALSE] %*% beta
   l <- crossprod(e)
   g <- t(-2*t(XTrain[index,,drop=FALSE]) %*% e)
   -g/length(index)
}

### create random data
set.seed(1)
N <- 1000
x <- rnorm(N)
X <- cbind(1, x)
y <- 1 + x + rnorm(N)
## training-validation
iTrain <- sample(N, 0.8*N)
XTrain <- X[iTrain,,drop=FALSE]
XValid <- X[-iTrain,,drop=FALSE]
yTrain <- y[iTrain]
yValid <- y[-iTrain]
cat("Analytic solution (training data):\n")
start <- c(const=10, x=10)
b0 <- drop(solve(crossprod(XTrain)) %*% crossprod(XTrain, yTrain))
names(b0) <- names(start)
tol <- 1e-1  # coefficient tolerance

## ---------- 1. working example
res <- maxSGA(loglik, gradlik, start=start,
            control=list(printLevel=0, iterlim=200,
                         SGA_batchSize=100, SGA_learningRate=0.1,
                         storeValues=TRUE),
            nObs=length(yTrain))
expect_equal(coef(res), b0, tolerance=tol)
                           # SGA usually ends with gradient not equal to 0 so we don't test that

## ---------- 2. full batch
res <- maxSGA(loglik, gradlik, start=start,
            control=list(printLevel=0, iterlim=200,
                         SGA_batchSize=NULL, SGA_learningRate=0.1,
                         storeValues=TRUE),
            nObs=length(yTrain))
expect_equal(coef(res), b0, tolerance=tol)

## ---------- 3. no gradient
res <- maxSGA(loglik, start=start,
              control=list(iterlim=1000, SGA_learningRate=0.02), nObs=length(yTrain))
expect_equal(coef(res), b0, tolerance=tol)
