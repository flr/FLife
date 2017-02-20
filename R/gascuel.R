#' @title Gascuel growth curve
#' 
#' @description 
#' Gascuel growth equation
#' 
#' @param age FLQuant, FLPar or numeric with ages 
#' @param params \code{FLPar}
#' @param ... any other arguments
#'  
#' @aliases gascuel gascuel-method gascuel,FLQuant,missing-method gascuel,FLPar,FLPar-method gascuel,FLPar,FLQuant-method
#' gascuel,FLPar,missing-method 
#' gascuel,numeric,FLQuant-method 
#' gascuel,numeric,numeric-method 
#' gascuel,numeric,FLPar-method 
#' gascuel,numeric,missing-method
#' gascuel,missing,FLPar-method 
#' gascuel,missing,missing-method
#' 
#' @return Returns a class of same type as \code{age} e.g. \code{FLQuant}

#' @export
#' @docType methods
#' @rdname gascuel
#' 
#' @seealso \code{\link{gompertz}}, \code{\link{vonB}} , \code{\link{richards}}
#'
#' @details 
#' Gascuel D., Fonteneau, A., and Capisano, C. (1992).
#' Modelisation d'une croissance en deux stances chez  #l'albacore (Thunnus albacares) de l'Atlantique Est. 
#' Aquat. Living Resour. 5: 155-172.
#'  
#' @examples
#' \dontrun{
#' gascuel(10)
#' }
#
setMethod("gascuel", signature(age="FLQuant",params="FLPar"),
          function(age,params,...){   
            res=gascuelFn(age,params)
            #res@units=""
            res})
setMethod("gascuel", signature(age="FLPar",params="FLPar"),
          function(age,params,...){   
            res=gascuelFn(age,params)
            #res@units=""
            res})
setMethod("gascuel", signature(age="numeric",params="FLPar"),
          function(age,params,...) 
            gascuelFn(age,params))
setMethod("gascuel", signature(age="FLQuant",params="missing"),
          function(age,params,...) {
            gascuelFn(age,params=FLPar(a=37.8,b=8.93,c=137.0,
                                       d=8.93,e=0.808,f=7.49))})
setMethod("gascuel", signature(age="FLPar",params="missing"),
          function(age,params,...) {
            gascuelFn(age,params=FLPar(a=37.8,b=8.93,c=137.0,
                                       d=8.93,e=0.808,f=7.49))})
setMethod("gascuel", signature(age="numeric",params="missing"),
          function(age,params,...) {
            gascuelFn(age,params=FLPar(a=37.8,b=8.93,c=137.0,
                                       d=8.93,e=0.808,f=7.49))})
setMethod("gascuel", signature(age="missing",params="FLPar"),
          function(age,params,length,...){   
            res=invGascuelFn(length,params,...)
            res})
setMethod("gascuel", signature(age="missing",params="missing"),
          function(age,params,length,...){   
            res=invGascuelFn(length,params=FLPar(a=37.8,b=8.93,c=137.0,
                                                 d=8.93,e=0.808,f=7.49),...)
            res})


# .expr2 <- 1-exp(t)
# .value <- .expr2^f
# 
# -(1-exp(t)^(f-1)*(f *exp(t)))

gascuelFn=function(age,params)
  params["a"]+
  params["b"]*age+
  (params["c"]-params["d"]*age)*
  (1-exp(- params["e"]*age))^params["f"]

deriv(l~a+b*age+(c-d*age)*(1-exp(-e*age))^f,"age")

dldt<-function(age,params){
  c=params["c"]
  b=params["b"]
  d=params["d"]
  e=params["e"]
  f=params["f"]
  
  .expr4 <- c - d * age
  .expr7 <- exp(-e * age)
  .expr8 <- 1 - .expr7
  .expr9 <- .expr8^f
  
  b+(.expr4*(.expr8^(f-1)*(f*(.expr7*e)))-d*.expr9)
}

invGascuelFn<-function(length,
                       params=FLPar(a=37.8,b=8.93,c=137.0,
                                    
                                    d=8.93,e=0.808,f=7.49),
                       age_limits=c(0,15),
                       timing    =0.5,
                       tol       =0.000001) {  
  
  names(params)=tolower(names(params))
  
  fn<-function(age,length,params,timing)
    (length-gascuel(age-timing,params))^2
  
  age=aaply(length, 1, function(x) 
    optimize(fn, age_limits,length=x,params=params,timing=-timing)$minimum)
  
  names(age)=names(length)
  
  age}

sliceGascuel<-function(length,
                       params=FLPar(a=37.8,b=8.93,c=137.0,
                                    d=8.93,e=0.808,f=7.49), 
                       age_limits=c(0,15),
                       timing=0.5,
                       tol   =0.000001) {  
  
  age=invGascuelFn(length,params,age_limits,timing,tol) 
  
  age=pmax(pmin(age, age_limits[2]), age_limits[1])
  age=as.integer(age)
  age}

