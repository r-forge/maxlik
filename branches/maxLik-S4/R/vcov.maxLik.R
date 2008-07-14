## maxLik
vcov.maxLik <- function(object, eigentol=1e-9, ...) {
   ## otherwise invert hessian
   activePar <- activePar(object)
   if(min(abs(eigen(hessian(object)[activePar,activePar],
                    symmetric=TRUE, only.values=TRUE)$values)) > eigentol) {
      varcovar <- matrix(0, nParam(object), nParam(object))
      varcovar[activePar,activePar] <-
          solve(-hessian(object)[activePar,activePar])
   }
   else
       varcovar <- NULL
   varcovar
}
