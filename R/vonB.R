#' @title von Bertalanffy growth curve
#'
#' @description 
#' Von Bertalanffy growth equation
#' 
#' @param age FLQuant, FLPar or numeric object with values corresponding to ages 
#' @param params \code{FLPar} object with parameters \code{linf}, \code{k} and \code{t0}
#' @param ... other arguments
#' 
#' @aliases vonB vonB-method vonB,FLPar,FLPar-method vonB,FLQuant,FLPar-method vonB,FLQuant,numeric-method vonB,missing,FLPar-method vonB,numeric,numeric-method vonB,numeric,FLPar-method
#' 
#' @return Returns an object of same class as \code{age}  e.g. \code{FLQuant}
#' 
#' @export
#' @docType methods
#' @rdname vonB
#' 
#' 
#' @seealso \code{\link{gompertz}}, \code{\link{gascuel}}, \code{\link{richards}}
#' 
#' @examples
#' \dontrun{
#' params=FLPar(linf=100,t0=0,k=.4)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' len=vonB(age,params)
#' 
#' #inverse growth curve
#' vonB(params=params,length=len)
#' }
setMethod("vonB", signature(age="FLQuant",params="FLPar"),
          function(age,params,...){   
            res=vonBFn(age,params)
            res@units=""
            res})
setMethod("vonB", signature(age="FLPar",params="FLPar"),
          function(age,params,...){   
            res=vonBFn(age,params)
            res@units=""
            res})
setMethod("vonB", signature(age="numeric",params="numeric"),
          function(age,params,...) 
            vonBFn(age,params))
setMethod("vonB", signature(age="FLQuant",params="numeric"),
          function(age,params,...) { 
            res=vonBFn(FLPar(params),age)
            res@units=""
            res})
setMethod("vonB", signature(age="numeric",params="numeric"),
          function(age,params,...) { 
            res=c(vonB(FLQuant(age),FLPar(params)))
            res})
setMethod("vonB", signature(age="numeric",params="FLPar"),
          function(age,params,...) { 
            res=c(vonB(FLQuant(age),params))
            res})

setMethod("vonB", signature(age="missing",params="FLPar"),
          function(age,params,length,...){  
            res=invVonBFn(length=length,params=params)
            res@units=""
            res})

vonBFn=function(age,params){
  
  dimnames(params)[[1]]=tolower(dimnames(params)[[1]])
  
  res=params["linf"]%*%(1.0-exp((-params["k"])%*%(age%-%params["t0"])))
  
  dimnames(res)[1:length(dim(res))]=dimnames(age)[1:length(dim(res))]
  res}

invVonBFn=function(length,params){
  res=log(1-(length%/%params["linf"]))%/% (-params["k"])%+%params["t0"]
  
  dimnames(res)=dimnames(length)
  res}

# library(numDeriv)
# params=FLPar(linf=318.9,k=0.093,t0=-0.970)
  # fnL=function(len) invVonB(params,FLQuant(len))
  # fnA=function(age)    vonB(params,FLQuant(age,dimnames=list(age=age)))
  # 
  # grad(fnL,fnA(15))

#vonB(FLPar(linf=100,k=.3,t0=-0.1),FLQuant(1))

