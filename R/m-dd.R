#' mdd
#' 
#'
#' Lorenzen natural mortality relationship where M is a function of weight, modified to 
#' explicitly included M as a function of numbers in a cohort, i.e. density dependence
#'
#'  @details
#'    
#' The Lorenzen natural mortality relationship is a function of mass-at-age i.e. 
#'    M=a*wt^b 
#' 
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
#'    M=aL/(1+exp(-k(n-ref)))*wt^b;
#'    
#' @param wt  mass at which M is to be predicted
#' @param params an \code{FLPar} with two values; i.e. a equal to M at unit mass 
#' and b a power term; defaults are a=0.3 and b=-0.288
#' @param scale, e.g. stock numbers now relative to a reference level, e.g. at virgin biomass. 
#' @param k steepness of relationship
#' 
#' @aliases mdd,FLQuant,FLPar-method
#' 
#' @export
#' @docType methods
#' @rdname mdd
#' 
#' @seealso \code{\link{lorenzen},\link{matdd}}
#'  
#' @examples
#' \dontrun{
#' library(FLBRP)
#' library(FLife)
#' 
#' data(ple4)
#' 
#' eql=brp(FLBRP(ple4))
#' fbar(eql)=FLQuant(c(seq(0,4,length.out=101)*refpts(eql)["msy","harvest"]))
#' stk=as(eql,"FLStock")
#' 
#' fbar(eql)=fbar(eql)[,1] 
#' scale=(stock.n(stk)%-%stock.n(eql))%/%stock.n(eql)
#' par=FLPar(m1=.2,m2=-0.288)
#' 
#' m=mdd(stock.wt(stk),par,scale)
#' ggplot(as.data.frame(m,drop=TRUE))+geom_line(aes(age,data,group=year,col=factor(year)))+theme(legend.position="none")
#' }
setGeneric('mdd', function(wt,params,...)
  standardGeneric('mdd'))

mddFn=function(wt,params,scale,k=1){
  
  map=(2/(1+exp(-k*scale)))
  
  (map%*%params["m1"])%*%(wt%^%params["m2"])}

setMethod('mdd', signature(wt='FLQuant',params='FLPar'),
          function(wt,params,scale,k=1) { 
            res=mddFn(wt,params,scale,k)
            res})