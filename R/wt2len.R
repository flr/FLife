#' @title Mass to length conversion
#'
#' @description 
#' converts weight to length
#' 
#' @param wt FLQuant, FLPar or numeric with length 
#' @param params \code{FLPar} with a and b, i.e. condition and scaling factors 
#' @param ... any other arguments
#' 
#' @aliases wt2len wt2len-method wt2len,FLCohort,FLPar-method wt2len,FLQuant,FLPar-method wt2len,numeric,FLPar-method
#' 
#' @return Returns an object of same class as \code{wt}  e.g. \code{FLQuant}
#' 
#' @export
#' @docType methods
#' @rdname wt2len
#' 
#' @seealso \code{\link{len2wt}}  
#' 
#' @examples
#' \dontrun{
#' params=FLPar(a=0.1,b=3)
#' wt2len(FLQuant(10),params)
#' }
#' 

#wt=a*exp(log(len)*b)
#a*exp(log(len)*b)=wt
#exp(log(len)*b)=(wt/a)
#log(len)*b=log(wt/a)
#log(len)=log(wt/a)/b
#len=exp(log(wt/a)/b)

setMethod("wt2len", signature(wt="FLQuant", params="FLPar"),function(wt,params,...) 
    exp(log(wt%/%params["a"])/params["b"]))
setMethod("wt2len", signature(wt="FLCohort",params="FLPar"),function(wt,params,...)  
    exp(log(wt%/%params["a"])/params["b"]))
setMethod("wt2len", signature(wt="numeric", params="FLPar"), function(wt,params,...) 
    exp(log(wt%/%params["a"])/params["b"]))