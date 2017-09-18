#' @title Random noise with different frequencies
#' 
#' @description A noise generator
#' 
#' @param n number of iterations
#' @param len an \code{FLQuant}
#' @param sd standard error for simulated series
#' @param b autocorrelation parameter a real number in [0,1] 
#' @param burn gets rid of 1st values i series
#' @param trunc get rid of values > abs(trunc)
#' @param what returns time series for year, cohort or age"
#' @param ... anyl
#' @aliases rnoise rnoise-method rnoise,numeric,FLQuant-method rnoise,numeric,missing-method
#' @aliases rlnoise rlnoise-method rlnoise,numeric,FLQuant-method rlnoise,numeric,missing-method
#' 
#' 
#' @export
#' @docType methods
#' @rdname rnoise
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
#' red <- rnoise(1000,flq,sd=.3,b=0.7)
#' plot(red)
#' acf(red)
#' 
#' data(ple4)
#' res=rnoise(1000,flq,sd=.3,b=0)
#' 
#' ggplot()+
#' geom_point(aes(year,age,size= data),
#'             data=subset(as.data.frame(res),data>0))+
#' geom_point(aes(year,age,size=-data),
#'             data=subset(as.data.frame(res),data<=0),colour="red")+
#' scale_size_area(max_size=4, guide="none")+
#' facet_wrap(~iter)
#' 
#' res=rnoise(4,m(ple4),burn=10,b=0.9,cohort=TRUE)
#' ggplot()+
#' geom_point(aes(year,age,size= data),
#'           data=subset(as.data.frame(res),data>0))+
#' geom_point(aes(year,age,size=-data),
#'           data=subset(as.data.frame(res),data<=0),colour="red")+
#' scale_size_area(max_size=4, guide="none")+
#' facet_wrap(~iter)
#' 
#' }
# n  =100
# len=FLQuant(0,dimnames=list(year=1:55)) 
# b  =0
# sd =0.3
# trunc=0
setMethod("rnoise", signature(n='numeric', len="FLQuant"),
    function(n=n,len=len,sd=0.3,b=0,burn=0,trunc=0,what=c("year","cohort","age")) {
      len=propagate(len,n)
      switch(what[1],
             "cohort"={object=as(len,"FLCohort")
                      res   =apply(object,c(2:6), function(x) 
                           t(noiseFn(length(x),sd,b,burn,trunc)))
                      res   =array(res,unlist(laply(dimnames(object),length)),
                            dimnames=dimnames(object))
                      res   =as(FLCohort(res),"FLQuant")
                      },
             "year"  ={res=apply(len,c(1,3:6), function(x) noiseFn(length(x),sd,b,burn,trunc))
                       res=as.FLQuant(res,dimnames=dimnames(len))},
             "age"   ={res=apply(len,c(2:6), function(x) noiseFn(length(x),sd,b,burn,trunc))
                       res=as.FLQuant(res,dimnames=dimnames(len))}
             )
      
      len+res})

setMethod("rnoise", signature(n='numeric', len="missing"),
          function(n=n,len=len,sd=0.3,b=0,burn=0,trunc=0,what=c("year","cohort","age")) {
             noiseFn(n,sd,b,burn,trunc)})

setMethod("rlnoise", signature(n='numeric', len="FLQuant"),
        function(n=n,len=len,sd=0.3,b=0,burn=0,trunc=0,what=c("year","cohort","age")) {
          exp(rnoise(n,len,sd,b,burn,trunc,what))})
                      
noiseFn<-function(len,sd=0.3,b=0,burn=0,trunc=0){
  
  if (burn<0) error("burn must be >=0")
  burn=burn+1
  x <- rep(0, len+burn) # going to hack off the first values at the end
  s <- rnorm(len+burn,mean=0,sd=sd)
  
  for(i in (1:(len+burn-1))){
    x[i+1] <- b * x[i] + s[i] * sqrt(1 - b^2)
    if(trunc>0){
      if (x[i+1] > (1-trunc))  x[i+1] <- ( 1-trunc)
      if (x[i+1] < (-1+trunc)) x[i+1] <- (-1+trunc)}
  }
  
  if (burn<=0) return(x)
  
  x<-x[-(seq(burn))]
  
  return(x)}

# setMethod("noise", signature(n='numeric', len="missing"),
#           function(n,len))

if (FALSE){
  library(FLCore)
  library(ggplotFL)
  data(ple4)
  res=noise(4,m(ple4),burn=10,b=0.9)
  ggplot()+
    geom_point(aes(year,age,size= data),
               data=subset(as.data.frame(res),data>0))+
    geom_point(aes(year,age,size=-data),
               data=subset(as.data.frame(res),data<=0),colour="red")+
    scale_size_area(max_size=4, guide="none")+
    facet_wrap(~iter)
}


## cohort effects
coEff=function(x,sd=1,b=0){
  
  dev        =log(rlnorm(length(dmns$cohort),0,cv))   
  for(i in 2:(length(dev)))
    dev[i]=dev[i]+dev[i-1]*rho
  parC[var]=parC[var]*exp(dev)
  
  tmp=len2wt(parC,vonB(ages(FLCohort(m(stk))),parC[c("linf","t0","k")]))
  res=window(as(tmp,"FLQuant"),start=1,end=dims(m(stk))$year)
  
  res}


cEff=function(stk,par,sigma,rho=0,var="k"){
  
  cv=sigma #getS(sigma,rho)
  dmns       =dimnames(par)
  dmns$cohort=dimnames(FLCohort(m(stk)))$cohort
  parC       =FLPar(rep(c(par),length(dmns$cohort)),dimnames=dmns[c(1,3,2)])
  units(parC)=""  
  
  dev        =log(rlnorm(length(dmns$cohort),0,cv))   
  for(i in 2:(length(dev)))
       dev[i]=dev[i]+dev[i-1]*rho
  parC[var]=parC[var]*exp(dev)
  
  tmp=len2wt(parC,vonB(ages(FLCohort(m(stk))),parC[c("linf","t0","k")]))
  res=window(as(tmp,"FLQuant"),start=1,end=dims(m(stk))$year)
  
  res}


#Spectral analysis function
spectra <- function(x,fs=1,norm = FALSE, pl = TRUE,omit=-(1:5))
{
  # Pad x with zeroes to make it's length a power of 2, i.e. length should be 2^something
  # This makes the fft faster
  oldx <- x # keep for later
  if(norm == TRUE) x <- x - mean(x)
  nfft <- (2^ceiling(log2(length(x))))
  x[(length(x)+1): nfft] <- 0
  fftx <- fft(x)
  # It's symmetrical so throw away second half. Only first 1 + nfft points are unique
  NoUnPo <- ceiling((nfft+1)/2) # Number of unique points
  fftx <- fftx[1:NoUnPo]
  # First element is DC component, last is the Nyquist component
  
  # Take magnitude of fft of x and scale the fft so that it is not a function of length
  mx <- abs(fftx) / length(x)
  # Take square of magnitude
  mx <- mx^2
  
  # As we dropped the first half of fft so multiply by 2 to keep same energy
  # DC component (first element) and Nyquist component (last element if even
  # number points (should be)) should not be multiplied by 2
  mx[2:(length(mx)-1)] <- mx[2:(length(mx)-1)] * 2
  
  f <- seq(from =0, to = NoUnPo-1) * fs/nfft # frequency axis for plot

  return(as.data.frame(list(mx = mx, f = f)))}

#ggplot(spectra(x))+geom_line(aes(f,mx))


