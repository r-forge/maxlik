### Does maxControl stuff behave?
### no need to test it on CRAN, hence private test
### 
### test for:
### 1. create maxControl object
### 2. SGA_batchSize NULL
### 3. negative batch size
### 4. more than 1 batch size

library(maxLik)
library(testthat)

maxControl(tol=1e-4, lambdatol=1e-5, qrtol=1e-6, qac="marquardt",
           marquardt_lambda0=0.1, marquardt_lambdaStep=3, marquardt_maxLambda=1e10,
           nm_alpha=2, nm_beta=1, nm_gamma=4,
           sann_temp=5, sann_tmax=100, sann_randomSeed=1,
           SGA_learningRate=0.5, SGA_batchSize=10, SGA_clip=1000, SGA_momentum=0.9,
           SG_patience=7, SG_patienceStep=10,
           iterlim=10, printLevel=3)

maxControl(SGA_batchSize=NULL)
try(maxControl(SGA_batchSize=-1))
try(maxControl(SGA_batchSize=2:3))

maxControl(SGA_clip=NULL)
try(maxControl(SGA_clip=-1))
try(maxControl(SGA_clip=2:3))
