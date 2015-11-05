#' vonB
#'
#' Von Bertalanffy growth equation
#' 
#' @param age FLQuant, FLPar or numeric with ages 
#' @param param
#' @param length FLQuant, FLPar or numeric with length, if supplied as a named paramreter
#' instead  of age then calculates ages.
#' 
#' @return Depends on the value of \code{data} 
#' 
#' #' @export
#' @docType methods
#' @rdname vonB
#' 
#' @seealso \code{\code{\link{gascuel}}}  
#' 
#' @examples
#' \dontrun{
#' param=FLPar(linf=100,t0=0,k=.4)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' len=vonB(age,param)
#' age=vonB(param,length=len)
#' }
setGeneric('vonB', function(age,param,...)
  standardGeneric('vonB'))

vonBFn=function(age,param){
 
  dimnames(param)[[1]]=tolower(dimnames(param)[[1]])
  
  res=param["linf"]%*%(1.0-exp((-param["k"])%*%(age%-%param["t0"])))
  
  dimnames(res)[1:5]=dimnames(age)[1:5]
  res}

invVonBFn=function(length,param){
  res=log(1-(length%/%param["linf"]))%/% (-param["k"])%+%param["t0"]

  #dimnames(res)=dimnames(length)
  res}

setMethod("vonB", signature(age="FLQuant",param="FLPar"),
          function(age,param,...){   
            res=FLife:::vonBFn(age,param)
            res@units=""
            res})
setMethod("vonB", signature(age="FLPar",param="FLPar"),
          function(age,param,...){   
            res=FLife:::vonBFn(age,param)
            res@units=""
            res})
setMethod("vonB", signature(age="numeric",param="numeric"),
          function(age,param,...) 
            vonBFn(age,param))
setMethod("vonB", signature(age="FLQuant",param="numeric"),
          function(age,param,...) { 
            res=FLife:::vonBFn(FLPar(param),age)
            res@units=""
            res})

setMethod("vonB", signature(age="missing",param="FLPar"),
          function(age,param,length,...){  
            res=FLife:::invVonBFn(length=length,param=param)
            res@units=""
            res})

# library(numDeriv)
# param=FLPar(linf=318.9,k=0.093,t0=-0.970)
  # fnL=function(len) invVonB(param,FLQuant(len))
  # fnA=function(age)    vonB(param,FLQuant(age,dimnames=list(age=age)))
  # 
  # grad(fnL,fnA(15))

#vonB(FLPar(linf=100,k=.3,t0=-0.1),FLQuant(1))

