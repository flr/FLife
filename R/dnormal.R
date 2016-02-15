#' dnormal
#' 
#'
#' Double normal ogive
#' 
#' @param age FLQuant or FLCohort 
#' @param params \code{FLPar} with parameters \code{a1} age at maximum, \code{sl} SD for lefthand limb and 
#' \code{sr} SD for righthand limb. 
#' @param ... any other arguments
#' 
#' @aliases dnormal-method dnormal,FLQuant,FLPar-method dnormal,FLPar,FLPar-method dnormal,numeric,numeric-method dnormal,FLQuant,numeric-method
#'                     
#' @return Depends on the value of \code{age} 
#'  
#' @export
#' @docType methods
#' @rdname ages
#' 
#' @seealso \code{\link{vonB}}  
#' 
#' @examples
#' \dontrun{
#' params=FLPar(a1=1,sl=3,sr=500)
#' wt2len(params,FLQuant(10))
#' }
setGeneric('dnormal', function(age,params,...)
  standardGeneric('dnormal'))

setMethod("dnormal", signature(age="FLQuant",params="FLPar"),
function(age,params,...){   
  res=dnormalFn(age,params)
  res@units=""
  res})
setMethod("dnormal", signature(age="FLPar",params="FLPar"),
          function(age,params,...){   
            res=dnormalFn(age,params)
            res@units=""
            res})
setMethod("dnormal", signature(age="numeric",params="numeric"),
          function(age,params,...) 
            dnormalFn(age,params))
setMethod("dnormal", signature(age="FLQuant",params="numeric"),
          function(age,params,...) { 
            res=dnormalFn(FLPar(params),age)
            res@units=""
            res})

dnormalFn_<-function(age,params){
  pow <-function(a,b) a^b
  func<- function(age,a1,sl,sr){
    if (age < a1)
      return(pow(2.0,-((age-a1)/sl*(age-a1)/sl)))
    else
      return(pow(2.0,-((age-a1)/sr*(age-a1)/sr)))}
  
  sapply(age,func,params["a1"],params["sl"],params["sr"])}


dnormalFn<-function(age,params){
  
  a1=FLQuant(1,dimnames=dimnames(age))%*%params["a1"]
  s =FLQuant(1,dimnames=dimnames(age))%*%params["sl"]
  sr=FLQuant(1,dimnames=dimnames(age))%*%params["sr"]
  
  if (dims(age)$iter==1 &  dims(a1)$iter>1)
    age=propagate(age,dims(a1)$iter)
  
  s[age>=a1]=sr[age>=a1]
  
  res=2.0^(-((age%-%a1)%/%s%*%(age%-%a1)%/%s))}
