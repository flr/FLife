#' @title Random noise with different frequencies
#' 
#' @description A noise generator for lognormal errors
#' 
#' @name rlnoise
#' 
#' @param n number of iterations
#' @param len an \code{FLQuant}
#' @param sd standard error for simulated series
#' @param b autocorrelation parameter a real number in [0,1] 
#' @param burn gets rid of 1st values i series
#' @param trunc get rid of values > abs(trunc)
#' @param what returns time series for year, cohort or age"
#' @param ... any other parameters
#' 
#' @aliases rlnoise rlnoise-method rlnoise,numeric,FLQuant-method rlnoise,numeric,missing-method
#' @export rlnoise
#' 
#' @docType methods
#' @rdname rlnoise
#'
#' @importFrom methods is
#'
#' @return A \code{FLQuant} with autocorrelation equal to B.
#' 
#' @references Ranta and Kaitala 2001 Proc. R. Soc.
#' vt = b * vt-1 + s * sqrt(1 - b^2)
#' s is normally distributed random variable with mean = 0
#' b is the autocorrelation parameter
#' @export
#' 
#' @examples
#' \dontrun{
#' flq=FLQuant(1:100)
#' white <- rnoise(1000,flq,sd=.3,b=0)
#' plot(white)
#' acf(white)
#' 
#' red <- rlnoise(1000,flq,sd=.3,b=0.7)
#' plot(red)
#' acf(red)
#' 
#' data(ple4)
#' res=rnoise(1000,log(flq),sd=.3,b=0)
#' 
#' ggplot()+
#' geom_point(aes(year,age,size= data),
#'             data=subset(as.data.frame(res),data>0))+
#' geom_point(aes(year,age,size=-data),
#'             data=subset(as.data.frame(res),data<=0),colour="red")+
#' scale_size_area(max_size=4, guide="none")+
#' facet_wrap(~iter)
#' 
#' res=rlnoise(4,log(m(ple4)),burn=10,b=0.9,cohort=TRUE)
#' ggplot()+
#' geom_point(aes(year,age,size= data),
#'           data=subset(as.data.frame(res),data>0))+
#' geom_point(aes(year,age,size=-data),
#'           data=subset(as.data.frame(res),data<=0),colour="red")+
#' scale_size_area(max_size=4, guide="none")+
#' facet_wrap(~iter)
#' 
#' }
