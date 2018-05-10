#' lorenzen
#'
#' Lorenzen natural mortality relationship estimate M as a function of weight. 
#' M=a*wt^b; 
#' 
#' @param wt  mass at which M is to be predicted
#' @param params an \code{FLPar} with two values; i.e. a equal to M at unit mass 
#' and b a power term; defaults are a=0.3 and b=-0.288
#' @param ... any other arguments
#' 
#' @aliases lorenzen lorenzen-method lorenzen,FLQuant,FLPar-method lorenzen,FLQuant,missing-method lorenzen,FLQuant,numeric-method  lorenzen,numeric,missing-method
#' 
#' @export
#' @docType methods
#' @rdname lorenzen
#' 
#' @seealso \code{\link{gislason}}
#'  
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' mass=FLQuant(c( 1.90, 4.23, 7.47,11.48,16.04,20.96,26.07,31.22,
#'                36.28,41.17,45.83,50.20,54.27,58.03,61.48,64.62),
#'              dimnames=list(age=1:16))
#' lorenzen(mass)
#' }

setMethod('lorenzen', signature(wt='FLQuant',params='FLPar'),
          function(wt,params,...){   
            res=params[1]%*%(wt%^%params[2])
            res@units='yr^-1'
            res})
setMethod('lorenzen', signature(wt='FLQuant',params='missing'),
      function(wt,m1=.3,m2=-0.288,...) { 
          res=lorenzenFn(wt,m1=m1,m2=m2)
          res@units='yr^-1'
          res})
setMethod('lorenzen', signature(wt='FLQuant',params='numeric'),
      function(wt,params,...) { 
          res=params[1]*wt^params[2]
          res@units='yr^-1'
          res})
setMethod('lorenzen', signature(wt='numeric',params='missing'),
          function(wt,m1=.3,m2=-0.288,...) { 
            res=lorenzenFn(wt,m1=m1,m2=m2)
            res})

lorenzenFn<-function(wt,m1=.3,m2=-0.288){
  if ("FLPar"%in%is(m2)) res=wt%^%m2  else res=wt^m2
  if ("FLPar"%in%is(m1)) res=m1%*%res else res=m1*res
  res}

m1<-function(m,wt){
  
  fn<-function(x,wt,ref) sum((lorenzen(wt,m1=x)-m)^2)
  
  optimize(fn, c(0, 100), tol=0.0000001,wt=wt,ref=m)$minimum}
  
