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
#' params=FLPar(linf=100,t0=0,k=.4)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' len=vonB(params,age)
#' age=vonB(params,length=len)
#' }
vonBFn=function(params,x){
 
  dimnames(par)[[1]]=tolower(dimnames(par)[[1]])
  
  res=par["linf"]%*%(1.0-exp((-par["k"])%*%(x%-%par["t0"])))
  
  dimnames(res)=dimnames(x)
  res}

invVonBFn=function(params,x){
  res=log(1-(x%/%par["linf"]))%/% (-par["k"])%+%par["t0"]

  dimnames(res)=dimnames(x)
  res}

setGeneric('vonB', function(params,x,...)
  standardGeneric('vonB'))

setMethod("vonB", signature(params="FLPar",x="FLQuant"),
          function(params,x,...){   
            res=vonBFn(params,x)
            res@units=""
            res})
setMethod("vonB", signature(params="FLPar",x="FLPar"),
          function(params,x,...){   
            res=vonBFn(params,x)
            res@units=""
            res})
setMethod("vonB", signature(params="numeric",x="numeric"),
          function(params,x,...) 
            vonBFn(params,x))
setMethod("vonB", signature(params="numeric",x="FLQuant"),
          function(params,x,...) { 
            res=vonBFn(FLPar(par),x)
            res@units=""
            res})

setMethod("vonB", signature(params="missing",x="FLPar"),
          function(params,x,length,...){   
            res=invVonBFn(params=x,x=length)
            res@units=""
            res})
# library(numDeriv)
# params=FLPar(linf=318.9,k=0.093,t0=-0.970)
  # fnL=function(len) invVonB(params,FLQuant(len))
  # fnA=function(age)    vonB(params,FLQuant(age,dimnames=list(age=age)))
  # 
  # grad(fnL,fnA(15))

#vonB(FLPar(linf=100,k=.3,t0=-0.1),FLQuant(1))

