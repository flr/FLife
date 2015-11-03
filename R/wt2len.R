#' wt2len
#'
#' converts weight to length
#' 
#' @param wt FLQuant, FLPar or numeric with length 
#' @param param FLPar with a and b, i.e. condition and scaling factors 
#' 
#' @return Depends on the value of \code{wt} 
#' 
#' @export
#' @docType methods
#' @rdname wt2len
#' 
#' @seealso \code{\link{len2wt}}  
#' 
#' @examples
#' \dontrun{
#' param=FLPar(a=1,b=3)
#' wt2len(FLQuant(10),param)
#' }
#' 
#' if (!isGeneric("rgamma"))
# gamma
setGeneric('wt2len', function(wt,param,...) standardGeneric('wt2len'))

setMethod("wt2len", signature(wt="FLQuant",param="FLPar"),
   function(wt,param,...) (wt/param["a"])^(1/param["b"]))
setMethod("wt2len", signature(wt="FLCohort",param="FLPar"),
   function(wt,param,...) (wt/param["a"])^(1/param["b"]))
setMethod("wt2len", signature(wt="numeric",param="FLPar"),
   function(wt,param,...) (wt/param["a"])^(1/param["b"]))