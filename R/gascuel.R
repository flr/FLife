globalVariables(c("invGascuelFn"))

gascuelFn=function(params,age)
  params["gas.a"]+
  params["gas.b"]*age+
  (params["gas.c"]-params["gas.d"]*age)*
  (1-exp(- params["gas.e"]*age))^params["gas.f"]

invGascuel<-function(length,
                     params, age_limits=c(0,15),
                     timing=0.5,
                     tol   =0.000001) {  
  
  names(params)=tolower(names(params))
  
  fn<-function(age,length,params,timing)
    (length-gascuel(params,age+timing))^2
  
  age=optimize(fn, age_limits,length=length,params=params,timing=-timing)$minimum
  
  age}

invGascuel<-function(length,
                     params, age_limits=c(0,15),
                     timing=0.5,
                     tol   =0.000001) {  
  
  age=invGascuel(length,params,age_limits,timing,tol) 
  
  age=pmax(pmin(age, age_limits[2]), age_limits[1])
  age=as.integer(age)
  age}

#' gascuel
#'
#' Gascuel growth equation
#' 
#' @param age FLQuant, FLPar or numeric with ages 
#' @param params \code{FLPar}
#' @param ... any other arguments
#' 
#' @aliases gascuel-method gascuel,FLPar,FLPar-method gascuel,FLPar,FLQuant-method  gascuel,FLPar,missing-method gascuel,numeric,FLQuant-method gascuel,numeric,numeric-method 
#' 
#' @return Depends on the value of \code{data} 
#' 
#' #' @export
#' @docType methods
#' @rdname gascuel
#' 
#' @seealso \code{\link{gompertz}}
#'
#' @details 
#' Gascuel D., Fonteneau, A., and Capisano, C. (1992).
#' Modelisation d'une croissance en deux stances chez  #l'albacore (Thunnus albacares) de l'Atlantique Est. 
#' Aquat. Living Resour. 5: 155-172.
#'  
#' @examples
#' \dontrun{
#' params=FLPar(gas.a=37.8,gas.b=8.93,gas.c=137.0,
#'              gas.d=8.93,gas.e=0.808,gas.f=7.49)
#' len=gascuel(params,age)
#' }
#
setGeneric('gascuel', function(params,age,...)
  standardGeneric('gascuel'))
setMethod("gascuel", signature(params="FLPar",age="FLQuant"),
          function(params,age,...){   
            res=gascuelFn(params,age)
            res@units=""
            res})
setMethod("gascuel", signature(params="FLPar",age="FLPar"),
          function(params,age,...){   
            res=gascuelFn(params,age)
            res@units=""
            res})
setMethod("gascuel", signature(params="numeric",age="numeric"),
          function(params,age,...) 
            gascuelFn(params,age))
setMethod("gascuel", signature(params="numeric",age="FLQuant"),
          function(params,age,...) { 
            res=gascuelFn(FLPar(params),age)
            res@units=""
            res})
setMethod("gascuel", signature(params="FLPar",age="missing"),
          function(params,age,length,...){   
            res=invGascuelFn(params,length)
            res@units=""
            res})

