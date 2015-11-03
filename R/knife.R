#' knife
#'
#' knife edge ogive
#' 
#' @param age FLQuant, FLPar or numeric with ages 
#' @param param
#' @param length FLQuant, FLPar or numeric with length, if supplied as a named paramreter
#' instead  of age then calculates ages.
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
#' param=FLPar(linf=100,t0=0,k=.4)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' len=knife(age,param)
#' age=knife(param,length=len)
#' }
setGeneric('knife', function(age,param,...)
  standardGeneric('knife'))

setMethod("knife", signature(age="FLQuant",param="FLPar"),
          function(age,param,...){   
            res=knifeFn(age,param)
            res@units=""
            res})
setMethod("knife", signature(age="FLPar",param="FLPar"),
          function(age,param,...){   
            res=knifeFn(age,param)
            res@units=""
            res})
setMethod("knife", signature(age="numeric",param="numeric"),
          function(age,param,...) 
            knifeFn(age,param))
setMethod("knife", signature(age="FLQuant",param="numeric"),
          function(age,param,...) { 
            res=knifeFn(FLPar(param),age)
            res@units=""
            res})

setMethod("knife", signature(age="missing",param="FLPar"),
          function(age,param,length,...){   
            res=invknifeFn(param,length)
            res@units=""
            res})

knifeFn<-function(age,par){
    res=age
    
    res[age> c(par["a1"])][]=0
    res[age<=c(par["a1"])][]=1
    
    return(res)}
