#' logistic
#' 
#'
#' logistic function
#' 
#' @param par
#' @param age
#' 
#' #' @export
#' @docType methods
#' @rdname logistic
#' 
#' @seealso \code{\code{\link{gompertz}}}
#' 
#' @examples
#' \dontrun{
#' par=lhPar(FLPar(linf=100))
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' mat=logistic(par,age)
#' }
setGeneric('logistic', function(age,params,...)
  standardGeneric('logistic'))

pow<-function(a,b) a^b
logisticFn<-function(age,params) { #x,a50,ato95,asym=1.0){  
  
  res =params["asym"]%/%(1.0+pow(19.0,(params["a50"]%-%age)%/%params["ato95"]))
  res[is.na(res)]=0
  asym=FLQuant(1,dimnames=dimnames(age))%*%params["asym"]
  grt =(params["a50"]%-%age)%/%params["ato95"] >  5
  lss =(params["a50"]%-%age)%/%params["ato95"] < -5
  
  res[grt]=0
  res[lss]=asym[lss]
  
  dmns          =dimnames(res)
  names(dmns)[1]="age"
  dimnames(res) =dmns
  
  return(res)}

# setMethod('logistic', signature(age='FLQuant',params='numeric'),
#           function(age,...) { 
#             res=logisticFn(age,FLPar(params))
#             res@units='proportion'
#             res})
setMethod('logistic', signature(age='FLQuant',params='FLPar'),
          function(age,params,...) { 
            res=logisticFn(age,params)
            res@units='proportion'
            res})