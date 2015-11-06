#' sigmoid
#'
#' sigmoid ogive
#' 
#' @param age FLQuant, FLPar or numeric with ages 
#' @param param
#' 
#' @return Depends on the type of \code{age} 
#' 
#' @export
#' @docType methods
#' @rdname sigmoid
#' 
#' @details
#' 
#' The sigmoid ogive is an S-shaped or sigmoid curve or sigmoidal functions, 
#' Verhulst hypothesizes that small populations increase geometrically, because the supply of resources exceeds demand. Then, as supply and demand balance, population growth is constant. Finally, as demand exceeds supply, population growth decreases at the same rate that it had increased. Verhulst describes this process with an equation that enables him to predict when a population will reach any given size (see Verhulst's Figure):
#' 
#' @seealso \code{\code{\link{dnormal},\link{knife}}}  
#' 
#' @examples
#' \dontrun{
#' param=FLPar(linf=100,t0=0,k=.4)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' len=sigmoid(age,param)
#' age=sigmoid(param,length=len)
#' }
setGeneric('sigmoid', function(age,param,...)
  standardGeneric('sigmoid'))

setMethod("sigmoid", signature(age="FLQuant",param="FLPar"),
          function(age,param,...){   
            res=sigmoidFn(age,param)
            res@units=""
            res})
setMethod("sigmoid", signature(age="FLPar",param="FLPar"),
          function(age,param,...){   
            res=sigmoidFn(age,param)
            res@units=""
            res})
setMethod("sigmoid", signature(age="numeric",param="numeric"),
          function(age,param,...) 
            sigmoidFn(age,param))
setMethod("sigmoid", signature(age="FLQuant",param="numeric"),
          function(age,param,...) { 
            res=sigmoidFn(FLPar(param),age)
            res@units=""
            res})

pow<-function(a,b) a^b
sigmoidFn<-function(age,param) { #x,a50,ato95,asym=1.0){  

  if (!("asym"%in%dimnames(param)$params))
    param=FLCore:::rbind(param,FLPar("asym" =1, iter=dims(param)$iter))
  if (!("ato95"%in%dimnames(param)$params))
    param=FLCore:::rbind(param,FLPar("ato95"=1, iter=dims(param)$iter))

  res<-param["asym"]%/%(1.0+pow(19.0,(param["a50"]%-%age)%/%param["ato95"]))
  asym=FLQuant(1,dimnames=dimnames(age))%*%param["asym"]
  res[(param["a50"]%-%age)%/%param["ato95"] >  5]<-0
  res[(param["a50"]%-%age)%/%param["ato95"] < -5]<-asym[(param["a50"]%-%age)%/%param["ato95"] < -5]

  dmns=dimnames(res)
  names(dmns)[1]="age"
  dimnames(res)=dmns

  return(res)}
