#' vonB
#'
#' Von Bertalanffy growth equation
#' 
#' @param par
#' @param age
#' 
#' #' @export
#' @docType methods
#' @rdname vonB
#' 
#' @seealso \code{\code{\link{gompertz}}}  
#' 
#' @examples
#' \dontrun{
#' par=FLPar(linf=100,t0=0,k=.4)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' len=vonB(age,par)
#' age=vonB(par,length=len)
#' }
vonBFn=function(par,x){
  res=par["linf"]%*%(1.0-exp((-par["k"])%*%(x%-%par["t0"])))
  
  dimnames(res)=dimnames(x)
  res}

invVonBFn=function(par,x){
  res=log(1-(x%/%par["linf"]))%/% (-par["k"])%+%par["t0"]

  dimnames(res)=dimnames(x)
  res}

setGeneric('vonB', function(par,x,...)
  standardGeneric('vonB'))

setMethod("vonB", signature(par="FLPar",x="FLQuant"),
          function(par,x,...){   
            res=vonBFn(par,x)
            res@units=""
            res})
setMethod("vonB", signature(par="FLPar",x="FLPar"),
          function(par,x,...){   
            res=vonBFn(par,x)
            res@units=""
            res})
setMethod("vonB", signature(par="numeric",x="numeric"),
          function(par,x,...) 
            vonBFn(par,x))
setMethod("vonB", signature(par="numeric",x="FLQuant"),
          function(par,x,...) { 
            res=vonBFn(FLPar(par),x)
            res@units=""
            res})
setMethod("vonB", signature(par="missing",x="FLPar"),
          function(par,x,length,...){   
            res=invVonBFn(par=x,x=length)
            res@units=""
            res})
# library(numDeriv)
# par=FLPar(linf=318.9,k=0.093,t0=-0.970)
  # fnL=function(len) invVonB(par,FLQuant(len))
  # fnA=function(age)    vonB(par,FLQuant(age,dimnames=list(age=age)))
  # 
  # grad(fnL,fnA(15))

#vonB(FLPar(linf=100,k=.3,t0=-0.1),FLQuant(1))

