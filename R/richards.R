#' @title Richards growth curve
#'
#' @description 
#' Richards growth equation
#' 
#' @param age FLQuant, FLPar or numeric object with values corresponding to ages 
#' @param params \code{FLPar} object with parameters \code{linf}, \code{k} and \code{t0}
#' @param ... other arguments
#' 
#' @aliases richards richards-method 
#'          richards,FLPar,FLPar-method 
#'          richards,FLQuant,FLPar-method 
#'          richards,FLQuant,numeric-method 
#'          richards,missing,FLPar-method 
#'          richards,numeric,numeric-method
#'          richards,numeric,FLPar-method    
#' 
#' @return Returns an object of same class as \code{age}  e.g. \code{FLQuant}
#' 
#' @export
#' @docType methods
#' @rdname richards
#' 
#' @seealso \code{\link{vonB}}, \code{\link{gompertz}}, \code{\link{gascuel}}
#' 
#' @examples
#' \dontrun{
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' len=richards(age,FLPar(linf=100,k=.4,b=.1,m=2))
#' }
setMethod("richards", signature(age="FLQuant",params="FLPar"),
          function(age,params,...){   
            res=richardsFn(age,params)
            res@units=""
            res})
setMethod("richards", signature(age="FLPar",params="FLPar"),
          function(age,params,...){   
            res=richardsFn(age,params)
            res@units=""
            res})
setMethod("richards", signature(age="numeric",params="FLPar"),
          function(age,params,...){   
            res=richardsFn(FLQuant(age,dimnames=list(age=seq(length(age)))),params)
            c(res)})
setMethod("richards", signature(age="numeric",params="numeric"),
          function(age,params,...) 
            richardsFn(age,params))
setMethod("richards", signature(age="FLQuant",params="numeric"),
          function(age,params,...) { 
            res=richardsFn(FLPar(params),age)
            res@units=""
            res})


# richards <- function(params,data) { #x, a50, ato95, sigma) {
#   beta <- ato95*log(19)/(log(2^sigma-1)-log((20/19)^sigma-1))
#   alpha <- a50+beta*log(2^sigma-1)/log(19)
#   
#   return((1/(1+19^(alpha-x)/beta))^1/sigma)} 


# bigeye FLPar(linf=178.63,k=0.424,b=-7.185,m=2880.4)
richardsFn=function(age,params){
  #linf/((1+exp(-kt+b))^m)
  
  dimnames(params)[[1]]=tolower(dimnames(params)[[1]])
  
  res=params["linf"]%/%exp(log(1+exp(-params["k"]%*%age+params["b"]))%*%params["m"])
  
  res}

invRichardsFn=function(length,params){
  res=log(1-(length%/%params["linf"]))%/% (-params["k"])%+%params["t0"]
  
  #dimnames(res)=dimnames(length)
  res}

# library(numDeriv)
# params=FLPar(linf=318.9,k=0.093,t0=-0.970)
  # fnL=function(len) invrichards(params,FLQuant(len))
  # fnA=function(age)    richards(params,FLQuant(age,dimnames=list(age=age)))
  # 
  # grad(fnL,fnA(15))

#richards(FLPar(linf=100,k=.3,t0=-0.1),FLQuant(1))

