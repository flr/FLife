#' dnormal
#'
#' Double normal ogive
#' 
#' @param params FLPar with parameters \code{a1} age at maximum, \code{sl} SD for lefthand limb and 
#' \code{sr} SD for righthand limb. 
#' @param age FLQuant or FLCohort 
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
#' params=FLPar(a1=1,sl=3,sr=500)
#' wt2len(params,FLQuant(10))
#' }
setGeneric('dnormal', function(params,age,...)
  standardGeneric('dnormal'))

setMethod("dnormal", signature(params="FLPar",age="FLQuant"),
function(params,age,...){   
  res=dnormalFn(params,age)
  res@units=""
  res})
setMethod("dnormal", signature(params="FLPar",age="FLPar"),
          function(params,age,...){   
            res=dnormalFn(params,age)
            res@units=""
            res})
setMethod("dnormal", signature(params="numeric",age="numeric"),
          function(params,age,...) 
            dnormalFn(params,age))
setMethod("dnormal", signature(params="numeric",age="FLQuant"),
          function(params,age,...) { 
            res=dnormalFn(FLPar(params),age)
            res@units=""
            res})

dnormalFn_<-function(params,age){
  pow <-function(a,b) a^b
  func<- function(age,a1,sl,sr){
    if (age < a1)
      return(pow(2.0,-((age-a1)/sl*(age-a1)/sl)))
    else
      return(pow(2.0,-((age-a1)/sr*(age-a1)/sr)))}
  
  sapply(age,func,params["a1"],params["sl"],params["sr"])}


dnormalFn<-function(params,age){
  
  a1=FLQuant(1,dimnames=dimnames(age))%*%params["a1"]
  s =FLQuant(1,dimnames=dimnames(age))%*%params["sl"]
  sr=FLQuant(1,dimnames=dimnames(age))%*%params["sr"]
  
  if (dims(age)$iter==1 &  dims(a1)$iter>1)
    age=propagate(age,dims(a1)$iter)
  
  s[age>=a1]=sr[age>=a1]
  
  res=2.0^(-((age%-%a1)%/%s%*%(age%-%a1)%/%s))}
