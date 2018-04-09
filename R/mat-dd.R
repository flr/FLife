#' matdd
#'
#' Logistic ogive for proportion mature-at-age, modified to explicitly included maturity 
#' as a function of numbers in a cohort, i.e. density dependence
#'
#' @details
#' The relationship can be explained by population density, since as fish grow they 
#' also die and so there is potentially less competition for resources between larger and 
#' older fish. Density dependence can be modelled by a logistic function, a sigmoid 
#' curve (or S shaped) curve, with equation 
#' 
#' f(x)=L/(1+exp(-k(x-x0)))
#'  
#' where
#' e  is the natural logarithm base (also known as Euler's number),
#' x0 is the x-value of the sigmoid's midpoint,
#' L  is the curve's maximum value, and
#' k  the steepness of the curve.
#' 
#' Combining the two functions gives
#'     
#'    O=aL/(1+exp(-k(n-ref)))*wt^b;
#'
#' @param age ages        
#' @param params an \code{FLPar} with two values; i.e. a equal to M at unit mass and b a power term; defaults are a=0.3 and b=-0.288
#' @param ... other arguments, such as scale, e.g. stock numbers now relative to a reference level, e.g. at virgin biomass and k steepness of relationship
#' 
#' @param scale reference 
#' @param k rate of change in density dependence
#' @param flagAge default is FALSE, i.e. density dependence is based on length ratherv than age
#' 
#' @aliases matdd matdd-method matdd,FLQuant,FLPar-method
#' 
#' @export
#' @docType methods
#' @rdname matdd
#' 
#' @seealso \code{\link{logistic},\link{mdd}}
#' 
#' @examples 
#' \dontrun{
#' library(FLRP)
#' library(FLife)
#' 
#' data(teleost)
#' par=teleost[,"Hucho hucho"]
#' par=lhPar(par)
#' hutchen=lhEql(par)
#' 
#' scale=stock.n(hutchen)[,25]%*%stock.wt(hutchen)
#' scale=(stock.n(hutchen)%*%stock.wt(hutchen)%-%scale)%/%scale
#' 
#' mat=matdd(ages(scale),par,scale,k=.5)   
#' 
#' ggplot(as.data.frame(mat))+
#'    geom_line(aes(age,data,col=factor(year)))+
#'    theme(legend.position="none")+
#'    scale_x_continuous(limits=c(0,15))
#' 
#'  }
setMethod('matdd', signature(age='FLQuant',params='FLPar'),
          function(age,params,scale,k=1,flagAge=TRUE) { 
            res=matddFn(age,params,scale,k,flagAge)
            res})

pow<-function(a,b) a^b

logisticFn2<-function(age,a50,ato95,asym) { #x,a50,ato95,asym=1.0){  
  
  res =asym%/%(1.0+pow(19.0,(a50%-%age)%/%ato95))
  dimnames(res)=dimnames(asym)
  res[is.na(res)]=0
  #asym=FLQuant(1,dimnames=dimnames(age))%*%asym
  
  for (i in dimnames(res)$iter){
    grt =iter(a50,i)%-%iter(age,i)%/%iter(ato95,i) >  5
    lss =iter(a50,i)%-%iter(age,i)%/%iter(ato95,i) < -5
    
    iter(res,i)[grt]=0
    iter(res,i)[lss]=iter(asym,i)[lss]
    }
  
  return(res)}

matddFn=function(age,params,scale,k=1,flagAge=FALSE){
  
  if (flagAge){
    a50 =FLPar(params["a50"]%*%(2/(1+exp(k*-scale))))
    asym=FLPar(FLQuant(c(params["asym"]),dimnames=dimnames(a50)))
  }else{
    asym=FLPar((2/(1+exp(k*scale))))    
    a50 =FLPar(FLQuant(c(params["a50"]),dimnames=dimnames(asym)))}
  
  ato95=(a50/a50)%*%params["ato95"]
    
  logisticFn2(age,a50,ato95,asym)}
