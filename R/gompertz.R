#' @title Gompertz growth equation
#' 
#' @description  
#' gompertz growth equation
#' 
#' @param age FLQuant, FLPar or numeric with ages 
#' @param params \code{FLPar} with parameters for \code{linf, a, k}
#' @param ... any other arguments
#' 
#' @aliases gompertz 
#' gompertz-method 
#' gompertz,FLQuant,FLPar-method g
#' gompertz,FLPar,FLPar-method gompertz,FLQuant,numeric-method gompertz,missing,FLPar-method gompertz,numeric,numeric-method
#' 
#' @return Returns an object of same class as \code{age} e.g. \code{FLQuant}
#' 
#' @exportMethod gompertz
#' @docType methods
#' @rdname gompertz
#'
#' @seealso \code{\link{gascuel}}, \code{\link{vonB}}, \code{\link{richards}}
#' 
#' @examples
#' \dontrun{
#' params=FLPar(linf=100,a=2,b=.4)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' gompertz(age,params)
#' }
#' 
setMethod("gompertz", signature(age="FLQuant",params="FLPar"),
          function(age,params,...){   
            res=gompertzFn(age,params)
            res@units=""
            res})
setMethod("gompertz", signature(age="FLPar",params="FLPar"),
          function(age,params,...){   
            res=gompertzFn(age,params)
            res@units=""
            res})
setMethod("gompertz", signature(age="numeric",params="numeric"),
          function(age,params,...) 
            gompertzFn(age,params))
setMethod("gompertz", signature(age="FLQuant",params="numeric"),
          function(age,params,...) { 
            res=gompertzFn(FLPar(params),age)
            res@units=""
            res})

setMethod("gompertz", signature(age="missing",params="FLPar"),
          function(age,params,length,...){   
            res=invgompertzFn(params,length)
            res@units=""
            res})

# bigeye FLPar(linf=179.13,k=0.4088,a=1.7268)                      

#linf*exp(-aexp(-kt))

gompertzFn=function(age,params) 
   params["linf"]%*%exp(-params["a"]%*%exp(log(params["k"])%*%age))


