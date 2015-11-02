#' cc
#'
#' Catch curve analysis
#' 
#' @param age
#' @param n
#' 
#' @return Depends on the value of \code{data} 
#' 
#' @export
#' @docType methods
#' @rdname cc
#' 
#' @seealso \code{\code{\link{powh}}}  
#' 
#' @examples
#' \dontrun{
#' params=FLPar(linf=100,t0=0,k=.4)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' len=vonB(params,age)
#' age=vonB(params,length=len)
#' }
setGeneric('cc', function(age,n,...)
  standardGeneric('cc'))

setMethod("cc", signature(age="numeric",n="numeric"),
          function(age,n,...){   
            res=ccFn(age,n)
            res})

setMethod("cc", signature(age="missing",n="FLQuant"),
          function(n){   
            dat=data.frame(n)
            res=with(dat,ccFn(age,data))
            res@units=""
            res})

setMethod("cc", signature(age="FLQuant",n="missing"),
          function(age){   
            dat=data.frame(age)
            res=with(dat,ccFn(age,data))
            res@units=""
            res})
ccFn=function(age,n){
  freq=freq/sum(freq)
  lm  =lm(log(freq)~age)
  hat =exp(predict(lm))
  sel =(freq/hat)/max(freq/hat)
  data.frame(age=age,obs=freq,hat=hat,sel=sel)}
