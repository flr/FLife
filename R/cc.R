#' cc
#'
#' Catch curve analysis
#' 
#' @param age age
#' @param n   frequency
#' @param ... any other arguments
#' 
#' @aliases cc cc-method cc,numeric,numeric-method cc,missing,FLQuant-method cc,FLQuant,missing-method
#' 
#' @return Depends on the value of \code{data} 
#' 
#' @export
#' @docType methods
#' @rdname cc
#' 
#' @seealso \code{\link{powh}}  
#' 
#' @examples
#' \dontrun{
#' params=FLPar(linf=100,t0=0,k=.4)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' len=vonB(params,age)
#' age=vonB(params,length=len)
#' }
setMethod("cc", signature(age="numeric",n="numeric"),
          function(age,n,...){  
            res=ccFn(age,n)
            res})

setMethod("cc", signature(age="missing",n="FLQuant"),
          function(age,n){   
            dat=data.frame(n)
            res=with(dat,ccFn(age,data))
            res@units=""
            res})

setMethod("cc", signature(age="FLQuant",n="missing"),
          function(age,n){   
            dat=data.frame(age)
            res=with(dat,ccFn(age,data))
            res@units=""
            res})
ccFn=function(age,n){
  lm  =lm(log(n)~age)
  hat =exp(predict(lm))
  sel =(n/hat)/max(n/hat)
  data.frame(age=age,obs=n,hat=hat,sel=sel)}
