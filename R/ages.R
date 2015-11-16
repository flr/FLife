#' ages
#'
#' Creates FLQuant/FLCohort with ages
#' 
#' @param object FLQuant or FLCohort 
#' 
#' @aliases ages-method ages,FLQuant-method ages,FLCohort-method
#' 
#' @return Depends on the value of \code{data} 
#'  
#' @export
#' @docType methods
#' @rdname ages
#' 
#' @seealso \code{\link{vonB}} 
#' 
#' @examples
#' \dontrun{
#' data(ple4)
#' ages(m(ple4))
#' }
setGeneric('ages', function(object, ...)
   standardGeneric('ages'))

setMethod("ages", signature(object="FLQuant"),
   function(object,timing=NULL){
      res<-FLQuant(dimnames(object)$age,dimnames=dimnames(object))

      if (is.null(timing))
         res<-sweep(res,4,(1:dim(res)[4]-1)/dim(res)[4],"+") else
         res<-sweep(res,4,timing,"+")

      return(res)})
setMethod("ages", signature(object="FLCohort"),
   function(object,timing=NULL){
      res<-FLCohort(dimnames(object)$age,dimnames=dimnames(object))

      if (is.null(timing))
         res<-sweep(res,4,(1:dim(res)[4]-1)/dim(res)[4],"+") else
         res<-sweep(res,4,timing,"+")

      return(res)})
