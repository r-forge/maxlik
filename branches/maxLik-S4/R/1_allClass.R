## All class definitions.  Read it first
library(methods)

setClass("lastStep",
         representation(theta0 = "numeric",
                        f0 = "numeric",
                        climb = "numeric"))

setClass("maxim",
         representation(maximum = "numeric",
                        estimate = "numeric",
                        gradient = "vector",
                        hessian = "matrix",
                        code = "integer",
                        message = "character",
                        iterations = "integer",
                        lastStep = "lastStep",
                        activePar = "logical",
                        type = "character"
                                        # a brief description of the maximization routine
                        ))

setClass("maxLik",
         representation("maxim"))

setClass("summary.maxLik",
         representation("maxLik",
                        results = "matrix",
                                        # adds a matrix of estimate, stds, t, Pr of coefficients
                        NActivePar = "integer")
                        )
