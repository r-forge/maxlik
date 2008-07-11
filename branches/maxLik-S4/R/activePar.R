## activePar: returns parameters which are free under maximisation (not fixed as constants)

activePar.maxim <- function(x)
    x@activePar
setMethod("activePar", "maxim", activePar.maxim)
rm(activePar.maxim)
