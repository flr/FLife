#' ages
#' 
#' Creates FLQuant/FLCohort with ages as entries
#' 
#' @param object \code{FLQuant} or \code{FLCohort} 
#' @param season the season dim
#' @param ... any other arguments
#' 
#' @aliases ages ages-method ages,FLQuant-method ages,FLCohort-method
#' 
#' @return \code{FLQuant} or \code{FLCohort} 
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
   function(object,season=NULL){
      res<-FLQuant(dimnames(object)$age,dimnames=dimnames(object))

      if (is.null(season))
         res<-sweep(res,4,(1:dim(res)[4]-1)/dim(res)[4],"+") else
         res<-sweep(res,4,season,"+")

      return(res)})
setMethod("ages", signature(object="FLCohort"),
   function(object,season=NULL){
      res<-FLCohort(dimnames(object)$age,dimnames=dimnames(object))

      if (is.null(season))
         res<-sweep(res,4,(1:dim(res)[4]-1)/dim(res)[4],"+") else
         res<-sweep(res,4,season,"+")

      return(res)})
