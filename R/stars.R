              
#' rod
#' 
#'
#' @description Regime shifts
#' Evidence for regime shifts are explored using a a sequential t-test algorithm 
#' (STARS; \cite{rodionov2004sequential}) as modified by Szuwalski et al., (submitted)
#' 
#' @param   object an object of class \code{FLQuant}
#' @param ... any other arguments
#' 
#' @aliases rod rod-method rod,FLQuant-method
#'
#' @export
#' @rdname rod
#'
#' @details Returns a data.frame 
#' 
#' @examples
#' \dontrun{
#'  #bug   
#'  object=rlnorm(1,FLQuant(0,dimnames=list(year=1:30)),.3)
#'  pg=rod(object) 
#'  plot(object) +
#'      geom_polygon(aes(year,data,group=regime),
#'      fill="lavender",col="blue",
#'      lwd=.25,data=pg,alpha=.2)
#' }
#' 
setMethod("rod", signature(object="FLQuant"),
          function(object, ...) 
            ddply(as.data.frame(object),.(iter), with, rodFn(data,year)))

#==Cody Szuwalski, 10/18,2012
#==this function identifies the first regime shift in a series given a significance level
#==based on a t-distribution (sig) and an assumed regime length (regLN) 
#==It must be used iteratively to identify all shifts in a time series (see example)

#==The key differences between this algorithm and STARS (Rodionov, 2004) is that 
#==when a new regime is detected, the algorithm begins searching at n+regimeLength,
#==where  is the first year of the regime identified and regimeLength is the 
#==assumed minimum length of a regime.  The original STARS begins searching at n+1.
#==Additionally, a 'Huber parameter' (which adjusts the influence of outliers)
#==is not included in this version of the algorithm.  

ROregimeSHFT=function(regLN,sig,series,shift=0){
   Shifts  =0
   startyr =1
   
   if(shift!=0) startyr=shift
 
   for (i in seq(startyr,(length(series)-(regLN)))){
     window=series[i:(regLN+(i-1))]
     basemean=mean(window,na.rm=T)
     basevar=var(window,na.rm=T)
     data=sig*sqrt(2*basevar/regLN)

   # the below line would be changed to something like:
   # if (series[i+1] > basemean+data | series[i+1] <basemean-data)
   # for the original algorithm (I have not dataed this)

  if (series[regLN+i] > basemean+data | series[regLN+i] <basemean-data){
    
     RSI=0
     if (series[regLN+i] > basemean+data){
       newmean=basemean+data
       RSI    =(series[regLN+i]-newmean)/(regLN*(sqrt(basevar)))
       down   =0
     }else{
       newmean=basemean-data
       RSI=(newmean-series[regLN+i])/(regLN*(sqrt(basevar)))
       down=1}
  
     origSIGN=sign(RSI)
     counter =min(regLN-1,(length(series)-(regLN+(i-1))))
     RSI     =0
   
     for (j in 1:counter){ 
       Shifts=0
       if(down==0)
            RSI= RSI + (series[regLN+i+j-1]-newmean)/(regLN*(sqrt(basevar)))
       else RSI= RSI + (newmean-series[regLN+i+j-1])/(regLN*(sqrt(basevar)))
    
       if (sign(RSI)!=origSIGN | RSI == Inf) break
                                              
       if (j == counter){
         Shifts=i
         break}
       }
   
    }
    
  if(Shifts!=0) break
  }
  
  return(Shifts)}

#==the 'Shifts' returned from "ROregimeSHFT" is the index (year) the length of 
#==the assumed regime before the regime occurs
#==that is, if Shift[1]==3, and the assumed regime length is 10
#==the regime shift occurs in year 13.

#==if ROregimeSHFT doesn't find a shift, it returns zero
#==the last entry in 'Shifts' (returned from the function) is either '0' or a value
#==for a year where there is a potential, but unconfirmed shift

#===============================================
#===example aplication=========================
#==============================================

#==an example time series==
# valMean=0
# valMean2=1
# valSD1=0.5
# valSD2=0.5
# val   =c(rnorm(12,valMean,valSD1),rnorm(12,valMean2,valSD2),rnorm(12,valMean,valSD1))			
# plot(val,type="l")
# 

rodFn=function(data,year=NULL,plot=FALSE){
  
  if (is.null(year)) year=seq(length(data))
  
  #==set the assumed minimum regime length==
  rgLN=10
  
  #== 'sig' should be based on an appropriate t-distribution given
  #== the number of observations in the series (but I didn't do that here)
  #== this significance level is roughly equal to p<0.1
  sig=1.68
  
  #==iteration for finding all the shifts in a time series==
  Shift=rep(NA,100)							# vector for recording shifts
   if(length(data)>rgLN)
   {
    try(Shift[1]<-ROregimeSHFT(regLN=rgLN,sig=1.68,series=data),TRUE)
    counts=2
    if(is.na(Shift[counts-1])==FALSE) 
    {
     while(Shift[counts-1]>0 & (Shift[counts-1]+rgLN+rgLN-1)<length(data))
      {
       Shift[counts]=ROregimeSHFT(rgLN,sig,data,(Shift[counts-1]+rgLN-1))
       counts=counts+1
      }
     }
    }   
  
  #===plot results=== 
  
  ShiftsVec=Shift[!is.na(Shift)]
  
  if (plot)
    plot(data,type="l")
  
  mn=NULL
  sd=NULL
  ln=NULL
  
  for(i in 1:length(ShiftsVec))
    {
    if(i==1){  
       regMean=mean(data[1:(Shift[i]+rgLN-1)])
       regSD=sd(data[1:(Shift[i]+rgLN-1)])
       mn=c(mn,regMean)
       sd=c(sd,regSD)
       ln=c(ln,Shift[i]+rgLN-1)
       
       if (plot)
         polygon(x=c(1,Shift[i]+rgLN-1,Shift[i]+rgLN-1,1),
       y=c(regMean+regSD,regMean+regSD,regMean-regSD,regMean-regSD),
       border=NA,col="#0000ff55")
      }
  
    if(i>1){
      endInd=ShiftsVec[i]+rgLN-1
      if(is.na(ShiftsVec[i+1]))
         {endInd=length(data)}
       
      regMean=mean(data[(Shift[i-1]+rgLN):endInd])
      regSD  =sd(  data[(Shift[i-1]+rgLN):endInd])
      mn=c(mn,regMean)
      sd=c(sd,regSD)
      ln=c(ln,endInd)
       
      if (plot)
        polygon(x=c(Shift[i-1]+rgLN,
                 endInd,
                 endInd,         
                 Shift[i-1]+rgLN),
             y=c(regMean+regSD,
                 regMean+regSD,
                 regMean-regSD,
                 regMean-regSD),
      border=NA,col="#0000ff55")}
      
  #print(ln)
  }

  ## regimes
  minyear=c(year[1],year[rev(rev(ln)[-1])]+1)
  maxyear=min(year)+ln-1
  
  dat=data.frame(mn=mn,sd=sd,i=factor(seq(length(mn))),
                 ln=ln,minyear=minyear,maxyear=maxyear)
  dat=data.frame(regime=dat$i,
                 data  =with(dat,c(mn+sd,  mn-sd,  mn-sd,  mn+sd)),
                 year  =with(dat,c(minyear,minyear,maxyear,maxyear)))
  #dat[do.call(order,dat),]
  }

