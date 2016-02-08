#' len2wt
#' 
#' 
#' converts length to weight
#' 
#' @param length age FLQuant, FLPar or numeric 
#' @param params \code{FLPar}
#' @param ... any other arguments
#' 
#' @aliases len2wt-method len2wt,FLCohort,FLPar-method len2wt,FLQuant,FLPar-method len2wt,numeric,FLPar-method
#' 
#' @return Depends on the value of \code{length} 
#' 
#' #' @export
#' @docType methods
#' @rdname len2wt
#' 
#' @seealso \code{\link{wt2len}} 
#' 
#' @examples
#' \dontrun{
#' params=FLPar(a=1,b=3)
#' len2wt(params,FLQuant(10))
#' }
## converts wt to len using condition factor
setGeneric('len2wt', function(length,params,...)
  standardGeneric('len2wt'))
  
setMethod("len2wt", signature(length="FLQuant",params="FLPar"),
   function(length,params) params["a"]*length^params["b"])
setMethod("len2wt", signature(length="FLCohort",params="FLPar"),
   function(length,params) params["a"]*length^params["b"])
setMethod("len2wt", signature(length="numeric",params="FLPar"),
   function(length,params) params["a"]*length^params["b"])
