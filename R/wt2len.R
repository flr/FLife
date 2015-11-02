#' wt2len
#'
#' converts weight to length
#' 
#' @param params FLPar with a and b, i.e. condition and scaling factors 
#' @param age FLQuant, FLPar or numeric with length 
#' 
#' @return Depends on the value of \code{data} 
#' 
#' #' @export
#' @docType methods
#' @rdname wt2len
#' 
#' @seealso \code{\code{\link{len2wt}}}  
#' 
#' @examples
#' \dontrun{
#' params=FLPar(a=1,b=3)
#' wt2len(params,FLQuant(10))
#' }
setGeneric('wt2len', function(params,data, ...)
  standardGeneric('wt2len'))
setMethod("wt2len", signature(params="FLPar", data="FLQuant"),
   function(params,data) (data/params["a"])^(1/params["b"]))
setMethod("wt2len", signature(params="FLPar", data="FLCohort"),
   function(params,data) (data/params["a"])^(1/params["b"]))
setMethod("wt2len", signature(params="FLPar", data="numeric"),
   function(params,data) (data/params["a"])^(1/params["b"]))