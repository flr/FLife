#' len2wt
#'
#' converts length to weight
#' 
#' @param length age FLQuant, FLPar or numeric 
#' @param param
#' 
#' @return Depends on the value of \code{length} 
#' 
#' #' @export
#' @docType methods
#' @rdname len2wt
#' 
#' @seealso \code{\code{\link{wt2len}}}  
#' 
#' @examples
#' \dontrun{
#' param=FLPar(a=1,b=3)
#' len2wt(param,FLQuant(10))
#' }
## converts wt to len using condition factor
setGeneric('len2wt', function(length,param,...)
  standardGeneric('len2wt'))
  
setMethod("len2wt", signature(length="FLQuant",param="FLPar"),
   function(length,param) param["a"]*length^param["b"])
setMethod("len2wt", signature(length="FLCohort",param="FLPar"),
   function(length,param) param["a"]*length^param["b"])
setMethod("len2wt", signature(length="numeric",param="FLPar"),
   function(length,param) param["a"]*length^param["b"])
