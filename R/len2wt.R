#' len2wt
#'
#' converts length to weight
#' 
#' @param params
#' @param length age FLQuant, FLPar or numeric 
#' 
#' @return Depends on the value of \code{data} 
#' 
#' #' @export
#' @docType methods
#' @rdname len2wt
#' 
#' @seealso \code{\code{\link{wt2len}}}  
#' 
#' @examples
#' \dontrun{
#' params=FLPar(a=1,b=3)
#' len2wt(params,FLQuant(10))
#' }
## converts wt to len using condition factor
setGeneric('len2wt', function(params,data,...)
  standardGeneric('len2wt'))
  
setMethod("len2wt", signature(params="FLPar", data="FLQuant"),
   function(params,data) params["a"]*data^params["b"])
setMethod("len2wt", signature(params="FLPar", data="FLCohort"),
   function(params,data) params["a"]*data^params["b"])
setMethod("len2wt", signature(params="FLPar", data="numeric"),
   function(params,data) params["a"]*data^params["b"])
