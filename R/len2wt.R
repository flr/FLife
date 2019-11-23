#' @title Length to weight conversion
#' 
#' @description  
#' Converts length to weight based on $W=aL^b$
#' 
#' @param length age FLQuant, FLPar or numeric 
#' @param params \code{FLPar}
#' @param ... any other arguments
#' 
#' @aliases len2wt len2wt-method len2wt,FLCohort,FLPar-method len2wt,FLQuant,FLPar-method len2wt,numeric,FLPar-method
#' 
#' @return Returns a class of same type as \code{length} e.g. \code{FLQuant}
#' @export
#' @docType methods
#' @rdname len2wt
#' 
#' @seealso \code{\link{wt2len}} 
#' 
#' @examples
#' \dontrun{
#' params=FLPar(a=1,b=3)
#' len2wt(FLQuant(10),params)
#' }
## converts wt to len using condition factor
setMethod("len2wt", signature(length="FLQuant",params="FLPar"),
   function(length,params) params["a"]%*%exp(log(length)%*%params["b"]))
setMethod("len2wt", signature(length="FLCohort",params="FLPar"),
   function(length,params) params["a"]%*%exp(log(length)%*%params["b"]))
setMethod("len2wt", signature(length="numeric",params="FLPar"),
   function(length,params) c(params["a"])*exp(log(length)*c(params["b"])))
