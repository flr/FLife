#' moment
#' 
#' @description 
#'   aaa 
#'          
#' @param object a vector holding a time series
#' @param n number of observations equal to length og object by default
#' @param na.rm \code{boolean} whether to remove NAs, TRUE by default 
#' @param ... any other arguments, i.e. x,n=rep(1,length(x)),na.rm=T
#' 
#' @aliases moment moment-method moment,numeric-method
#' @return a \code{vector} with the inter-annual variation each time step
#' @export
#' @docType methods
#' @rdname moment
#' 
#' @aliases moment-method moment,numeric-method 
#' 
#' @examples
#' x=1

setMethod("moment", signature(object='numeric'),
          function(object,n=rep(1,length(object)),na.rm=T) { 
  if(length(n)==1) n=rep(n,length(object)) 
  
  mn= sum(object*n,            na.rm=na.rm)/sum(n,na.rm=na.rm)
  sd=(sum(n*(object-mn)^2,     na.rm=na.rm)/sum(n,na.rm=na.rm))^.5
  sk= sum(n*((object-mn)/sd)^3,na.rm=na.rm)/sum(n,na.rm=na.rm)
  ku= sum(n*((object-mn)/sd)^4,na.rm=na.rm)/sum(n,na.rm=na.rm)-3
  
  ## weighted median
  n=unlist(c(n))
  object=unlist(c(object))
  
  cumn=cumsum(n)/sum(n)
  max.=max((1:length(n))[cumn<.50])
  min.=min((1:length(n))[cumn>=.50])
  
  me=(x[min.]*n[min.]+x[max.]*n[max.])/(n[min.]+n[max.])
  
  return(c(mn=mn,sd=sd,sk=sk,ku=ku,me=me))})

moments<-function(x,n,p=1) (sum(x^p*n)/sum(n))^(1/p)

