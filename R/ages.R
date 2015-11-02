#' ages
#'
#' Creates FLQuant/FLCohort with ages
#' 
#' @param age FLQuant or FLCohort 
#' 
#' @return Depends on the value of \code{data} 
#' 
#' @return Depends on the value of \code{data} 
#' 
#' #' @export
#' @docType methods
#' @rdname ages
#' 
#' @seealso \code{\code{\link{vonB}}}  
#' 
#' @examples
#' \dontrun{
#' params=FLPar(a=1,b=3)
#' wt2len(params,FLQuant(10))
#' }
setGeneric('ages', function(age, ...)
   standardGeneric('ages'))

setMethod("ages", signature(age="FLQuant"),
   function(age,timing=NULL){
      res<-FLQuant(dimnames(age)$age,dimnames=dimnames(age))

      if (is.null(timing))
         res<-sweep(res,4,(1:dim(res)[4]-1)/dim(res)[4],"+") else
         res<-sweep(res,4,timing,"+")

      return(res)})
setMethod("ages", signature(age="FLCohort"),
   function(age,timing=NULL){
      res<-FLCohort(dimnames(age)$age,dimnames=dimnames(age))

      if (is.null(timing))
         res<-sweep(res,4,(1:dim(res)[4]-1)/dim(res)[4],"+") else
         res<-sweep(res,4,timing,"+")

      return(res)})
