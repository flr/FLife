#' dnormal
#'
#' Double normal ogive
#' 
#' @param age FLQuant or FLCohort 
#' @param param FLPar with parameters \code{a1} age at maximum, \code{sl} SD for lefthand limb and 
#' \code{sr} SD for righthand limb. 
#' 
#' @return Depends on the value of \code{age} 
#'  
#' #' @export
#' @docType methods
#' @rdname ages
#' 
#' @seealso \code{\code{\link{vonB}}}  
#' 
#' @examples
#' \dontrun{
#' param=FLPar(a1=1,sl=3,sr=500)
#' wt2len(param,FLQuant(10))
#' }
setGeneric('dnormal', function(age,param,...)
  standardGeneric('dnormal'))

setMethod("dnormal", signature(age="FLQuant",param="FLPar"),
function(age,param,...){   
  res=dnormalFn(age,param)
  res@units=""
  res})
setMethod("dnormal", signature(age="FLPar",param="FLPar"),
          function(age,param,...){   
            res=dnormalFn(age,param)
            res@units=""
            res})
setMethod("dnormal", signature(age="numeric",param="numeric"),
          function(age,param,...) 
            dnormalFn(age,param))
setMethod("dnormal", signature(age="FLQuant",param="numeric"),
          function(age,param,...) { 
            res=dnormalFn(FLPar(param),age)
            res@units=""
            res})

dnormalFn_<-function(age,param){
  pow <-function(a,b) a^b
  func<- function(age,a1,sl,sr){
    if (age < a1)
      return(pow(2.0,-((age-a1)/sl*(age-a1)/sl)))
    else
      return(pow(2.0,-((age-a1)/sr*(age-a1)/sr)))}
  
  sapply(age,func,param["a1"],param["sl"],param["sr"])}


dnormalFn<-function(age,param){
  
  a1=FLQuant(1,dimnames=dimnames(age))%*%param["a1"]
  s =FLQuant(1,dimnames=dimnames(age))%*%param["sl"]
  sr=FLQuant(1,dimnames=dimnames(age))%*%param["sr"]
  
  if (dims(age)$iter==1 &  dims(a1)$iter>1)
    age=propagate(age,dims(a1)$iter)
  
  s[age>=a1]=sr[age>=a1]
  
  res=2.0^(-((age%-%a1)%/%s%*%(age%-%a1)%/%s))}
