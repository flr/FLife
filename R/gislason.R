#' gislason
#'
#' gislason natural mortality relatoinship estimate M as a function of weight. 
#' M=a*length^b; 
#' 
#' @param length  mass at which M is to be predicted
#' @param par an FLPar with two values; i.e. a equal to M at unit mass 
#' and b a power term; defaults are a=0.3 and b=-0.288
#'
#' @import FLCore 
#' 
#' @export
#' @docType methods
#' @rdname gislason
#' 
#' @seealso \code{\link{gislason}}
#'  
#' @examples
#' \dontrun{
#' length=FLQuant(c( 1.90, 4.23, 7.47,11.48,16.04,20.96,26.07,31.22,
#'                36.28,41.17,45.83,50.20,54.27,58.03,61.48,64.62),
#'              dimnames=list(age=1:16))
#' gislason(mass)
#' }
setGeneric('gislason', function(length,param,...)
  standardGeneric('gislason'))

gislasonFn<-function(length,param) {
  
  # Natural mortality parameters from Model 2, Table 1 gislason 2010
  if (!all(c("m1","m2")%in%dimnames(param)$params))
    param=FLCore:::rbind(param,FLPar(m1= 0.55*(param["linf"]^1.44)%*%param["k"], iter=dims(param)$iter),
                               FLPar(m2=-1.61                                  , iter=dims(param)$iter))
  
  param["m1"]%*%(exp(log(length)%*%param["m2"]))}

setMethod('gislason', signature(length='FLQuant',param='missing'),
      function(length,...) { 
          res=gislasonFn(length,param)
          res@units='yr^-1'
          res})
setMethod('gislason', signature(length='FLQuant',param='numeric'),
      function(length,param,...) { 
          res=gislasonFn(length,param)
          res@units='yr^-1'
          res})
setMethod('gislason', signature(length='FLQuant',param='FLPar'),
      function(length,param,...){   
          res=gislasonFn(length,param)
          res@units='yr^-1'
          res})

