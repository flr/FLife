#' knife
#' 
#'
#' knife edge ogive
#' 
#' @param age FLQuant, FLPar or numeric with ages 
#' @param params \code{FLPar}
#' @param ... any other arguments
#' 
#' @return Depends on the value of \code{data} 
#' 
#' @docType methods
#' @rdname knife
#' 
#' @details
#' The knife ogive is an S-shaped or knife curve or knifeal functions, 
#' Verhulst hypothesizes that small populations increase geometrically, because the supply of resources exceeds demand. Then, as supply and demand balance, population growth is constant. Finally, as demand exceeds supply, population growth decreases at the same rate that it had increased. Verhulst describes this process with an equation that enables him to predict when a population will reach any given size (see Verhulst's Figure):
#' 
#' @seealso \code{\link{dnormal},\link{sigmoid}}  
#' 
#' @examples
#' \dontrun{
#' params=FLPar(linf=100,t0=0,k=.4)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' len=knife(age,params)
#' age=knife(params,length=len)
#' }
setGeneric('knife', function(age,params,...)
  standardGeneric('knife'))

setMethod("knife", signature(age="FLQuant",params="FLPar"),
          function(age,params,...){   
            res=knifeFn(age,params)
            res@units=""
            res})
setMethod("knife", signature(age="FLPar",params="FLPar"),
          function(age,params,...){   
            res=knifeFn(age,params)
            res@units=""
            res})
setMethod("knife", signature(age="numeric",params="numeric"),
          function(age,params,...) 
            knifeFn(age,params))
setMethod("knife", signature(age="FLQuant",params="numeric"),
          function(age,params,...) { 
            res=knifeFn(FLPar(params),age)
            res@units=""
            res})
knifeFn<-function(age,par){
    res=age
    
    res[age> c(par["a1"])][]=0
    res[age<=c(par["a1"])][]=1
    
    return(res)}
