#' @title Sigmoid ogive
#' 
#' @description 
#' sigmoid ogive
#' 
#' @param age FLQuant, FLPar or numeric with ages 
#' @param params \code{FLPar} with age at whick probability is 0.5 (\code{a50}) and offset to 0.95 (\code{ato95})
#' @param ... any other arguments
#' 
#' @aliases sigmoid sigmoid-method sigmoid,FLPar,FLPar-method sigmoid,FLQuant,FLPar-method sigmoid,FLQuant,numeric-method sigmoid,numeric,numeric-method
#' 
#' @return returns an object of same type as \code{age} e.g. \code{FLQuant}
#' 
#' @export
#' @docType methods
#' @rdname sigmoid
#' 
#' @seealso \code{\link{knife}}, \code{\link{dnormal}}, \code{\link{logistic}}
#' 
#' @details
#' 
#' The sigmoid ogive is an S-shaped or sigmoid curve or sigmoidal functions, 
#' Verhulst hypothesizes that small populations increase geometrically, because the supply of resources exceeds demand. Then, as supply and demand balance, population growth is constant. Finally, as demand exceeds supply, population growth decreases at the same rate that it had increased. Verhulst describes this process with an equation that enables him to predict when a population will reach any given size (see Verhulst's Figure):
#' 
#' @seealso \code{\link{dnormal},\link{knife}}
#' 
#' @examples
#' \dontrun{
#' params=FLPar(a50=4,ato95=1)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' sigmoid(age,params)
#' }
setMethod("sigmoid", signature(age="FLQuant",params="FLPar"),
          function(age,params,...){   
            res=sigmoidFn(age,params)
            res@units=""
            res})
setMethod("sigmoid", signature(age="FLPar",params="FLPar"),
          function(age,params,...){   
            res=sigmoidFn(age,params)
            res@units=""
            res})
setMethod("sigmoid", signature(age="numeric",params="numeric"),
          function(age,params,...) 
            sigmoidFn(age,params))
setMethod("sigmoid", signature(age="FLQuant",params="numeric"),
          function(age,params,...) { 
            res=sigmoidFn(FLPar(params),age)
            res@units=""
            res})

sigmoidFn<-function(age,params) { #x,a50,ato95,asym=1.0){  

  if (!("asym"%in%dimnames(params)$params))
    params=rbind(params,FLPar("asym" =1, iter=dims(params)$iter))
  if (!("ato95"%in%dimnames(params)$params))
    params=rbind(params,FLPar("ato95"=1, iter=dims(params)$iter))

  res<-params["asym"]%/%(1.0+pow(19.0,(params["a50"]%-%age)%/%params["ato95"]))
  asym=FLQuant(1,dimnames=dimnames(age))%*%params["asym"]
  res[(params["a50"]%-%age)%/%params["ato95"] >  5]<-0
  res[(params["a50"]%-%age)%/%params["ato95"] < -5]<-asym[(params["a50"]%-%age)%/%params["ato95"] < -5]

  dmns=dimnames(res)
  names(dmns)[1]="age"
  dimnames(res)=dmns

  return(res)}
