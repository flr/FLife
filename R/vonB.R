#' vonB
#'
#' Von Bertalanffy growth equation
#' 
#' @param params
#' @param age FLQuant, FLPar or numeric with ages 
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
#' params=FLPar(linf=100,t0=0,k=.4)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' len=vonB(params,age)
#' age=vonB(params,length=len)
#' }
setGeneric('vonB', function(params,age,...)
  standardGeneric('vonB'))

vonBFn=function(params,age){
 
  dimnames(params)[[1]]=tolower(dimnames(params)[[1]])
  
  res=params["linf"]%*%(1.0-exp((-params["k"])%*%(age%-%params["t0"])))
  
  dimnames(res)[1:5]=dimnames(age)[1:5]
  res}

invVonBFn=function(params,age){
  res=log(1-(age%/%params["linf"]))%/% (-params["k"])%+%params["t0"]

  dimnames(res)=dimnames(age)
  res}

setMethod("vonB", signature(params="FLPar",age="FLQuant"),
          function(params,age,...){   
            res=vonBFn(params,age)
            res@units=""
            res})
setMethod("vonB", signature(params="FLPar",age="FLPar"),
          function(params,age,...){   
            res=vonBFn(params,age)
            res@units=""
            res})
setMethod("vonB", signature(params="numeric",age="numeric"),
          function(params,age,...) 
            vonBFn(params,age))
setMethod("vonB", signature(params="numeric",age="FLQuant"),
          function(params,age,...) { 
            res=vonBFn(FLPar(params),age)
            res@units=""
            res})

setMethod("vonB", signature(params="FLPar",age="missing"),
          function(params,age,length,...){   
            res=invVonBFn(params,length)
            res@units=""
            res})

# library(numDeriv)
# params=FLPar(linf=318.9,k=0.093,t0=-0.970)
  # fnL=function(len) invVonB(params,FLQuant(len))
  # fnA=function(age)    vonB(params,FLQuant(age,dimnames=list(age=age)))
  # 
  # grad(fnL,fnA(15))

#vonB(FLPar(linf=100,k=.3,t0=-0.1),FLQuant(1))

