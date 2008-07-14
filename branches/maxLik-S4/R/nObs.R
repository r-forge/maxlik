## Return #of observations for models

nObs.lm <- function(x)
    nrow(x$qr$qr)

nObs.default <- function(x)
    x$param$nObs

