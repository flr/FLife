#' @title Fills an \code{FLQuant} with ages
#' 
#' @description  
#' Creates \code{FLQuant} and \code{FLCohort} with ages as entries
#' 
#' @param object \code{FLQuant} or \code{FLCohort} 
#' @param ... any other arguments
#' 
#' @aliases ages ages-method ages,FLQuant-method ages,FLCohort-method
#' 
#' @return \code{FLQuant} or \code{FLCohort} depending on what the first argument was
#'  
#' @export
#' @docType methods
#' @rdname ages
#' 
#' @seealso \code{\link{knife}} \code{\link{gascuel}} \code{\link{sigmoid}} \code{\link{gompertz}} \code{\link{vonB}} \code{\link{dnormal}} \code{\link{logistic}}
#' 
#' @examples
#' \dontrun{
#' data(ple4)
#' ages(m(ple4))}
setMethod("ages", signature(object="FLQuant"),
   function(object){
      res<-FLQuant(dimnames(object)$age,dimnames=dimnames(object))

      if (dim(object)[4]>1)
         res<-sweep(res,4,(1:dim(res)[4]-1)/dim(res)[4],"+")
        
      return(res)})
setMethod("ages", signature(object="FLCohort"),
   function(object){
      res<-FLCohort(dimnames(object)$age,dimnames=dimnames(object))

      if (dim(object)[4]>1)
         res<-sweep(res,4,(1:dim(res)[4]-1)/dim(res)[4],"+") 
         
      return(res)})
