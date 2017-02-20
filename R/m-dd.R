#' mdd
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
#' @param params an \code{FLPar} with two values; i.e. a equal to M at unit mass and b a power term; defaults are a=0.3 and b=-0.288
#' @param scale reference 
#' @param k rate of change in density dependence
#' @param ... other arguments, such as scale, e.g. stock numbers now relative to a reference level, e.g. at virgin biomass and k steepness of relationship
#' 
#' @aliases mdd mdd-method mdd,FLQuant,FLPar-method
#' 
#' @export
#' @docType methods
#' @rdname mdd
#' 
#' @seealso \code{\link{lorenzen}}
#'  
#' @examples
#' \dontrun{
#' library(FLBRP)
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
#' m=mdd(wt2len(stock.wt(hutchen),par),params=par,scale,k=.9) 
#'  
#' ggplot(as.data.frame(m))+
#'    geom_line(aes(age,data,col=factor(year)))+
#'    theme(legend.position="none")+
#'    scale_x_continuous(limits=c(0,15))
#' 
#' m=mdd(stock.wt(hutchen),params=FLPar(m1=3,m2=-0.288),scale,k=1.2,m=lorenzen)   
#' }
setMethod('mdd', signature(object='FLQuant',params='FLPar'),
          function(object,params,scale,k=1,m=gislason) { 
            
            map=(2/(1+exp(-k*scale)))
            
            res=map%*%m(object,params)
            res})

