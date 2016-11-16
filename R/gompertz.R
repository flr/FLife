globalVariables(c("invgompertzFn"))

#' @title{gompertz growth equation}
#' @description  
#' gompertz growth equation
#' 
#' @param age FLQuant, FLPar or numeric with ages 
#' @param params \code{FLPar} with parameters for \code{linf, a, b}
#' @param ... any other arguments
#' 
#' @aliases gompertz 
#' gompertz-method 
#' gompertz,FLQuant,FLPar-method g
#' gompertz,FLPar,FLPar-method gompertz,FLQuant,numeric-method gompertz,missing,FLPar-method gompertz,numeric,numeric-method
#' 
#' @return Depends on the value of \code{data} 
#' 
#' @exportMethod gompertz
#' @docType methods
#' @rdname gompertz
#' 
#' @examples
#' \dontrun{
#' par=FLPar(linf=100,t0=0,k=.4)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' len=gompertz(par,age)
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

gompertzFn=function(params,age) 
   params["linf"]%*%exp(-params["a"]%*%params["b"]%^%age)


