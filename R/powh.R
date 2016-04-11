globalVariables(c("ddply",".","x"))

#' powh
#' 
#' @description 
#' 
#' Estimates growth and mortality parameters from length frequency data.
#'  
#' @param len vector with length distribution
#' @param n vector with numbers in each length bin
#' @param ... any other argument, i.e. weights =TRUE, fromMode =FALSE
#' 
#' @return a \code{data.frame} \code{mn} (mean), \code{diff} (difference), 
#' \code{len} (length) and \code{n} (frequency)
#'  
#' @details
#'  
#'  Beverton and Holt (1956) developed a method to estimate population parameters 
#'  such as total mortality (Z) from length data i.e.
#' 
#' \deqn{Z=K\frac{L_{infinity}-\overline{L}}{\overline{L}-L'}}
#'                             
#' Powell (1979) then developed a method, extended by Wetherall et al. (1987), 
#' to estimate growth and mortality parameters. This assumes that the right hand tail 
#' of a length frequency distribution is determined by the asymptotic length L 
#' and the ratio between Z and the growth rate K.
#' 
#' The Beverton and Holt methods assumes good estimates for K and $L_{infinity}$, 
#' while the Powell-Wetherall method only requires an estimate
#' of K, since $L_{infinity}$ is estimated by the method as well as Z/K. These method 
#' therefore provide estimates for Z/K, if K is unknown and 
#' Z if K is known.  
#' 
#' As well as assuming that growth follows the von Bertalanffy growth function, 
#' it is also assumed that the population is in a steady state 
#' with constant exponential mortality, no changes in selection pattern of the 
#' fishery and constant recruitment.
#' 
#' In the Powell-Wetherall method L' can take any value between the 
#' smallest and largest sizes. Equation 1 then provides a series of estimates of
#' Z and since 
#' 
#' \deqn{
#' \overline{L}-L'=a+bL'
#' }
#' 
#' a and b can be estimated by a regression analysis where 
#' 
#' \deqn{b={-K}/{Z+K}}
#' \deqn{a=-bL_{infinity}}
#' 
#' Therefore plotting \deqn{\overline{L}-L} against $L'$ provides an estimate 
#' of $L_{infinity}$ and Z/K from
#' 
#' \deqn{L_{infinity}=-a/b}
#' \deqn{Z/K={-1-b}/{b}}
#' 
#' If K is known Z can also be esimated
#'
#' @aliases powh-method powh,numeric,numeric-method 
#' 
#' @references
#' 
#'  R. Beverton and S. Holt. Review of method for estimating mortality rates in 
#'  exploited fish populations, with special reference to sources of bias in 
#'  catch sampling. Rapports et Proces-Verbaux., 140(1): 67--83, 1956.   
#'  
#'  D. G. Powell.
#'  Estimation of mortality and growth parameters from the length
#'  frequency of a catch [model].
#'  \emph{Rapports et Proces-Verbaux des Reunions}, 175, 1979.
#'  
#'  J. Wetherall, J. Polovina, and S. Ralston.
#'  Estimating growth and mortality in steady-state fish stocks from
#'  length-frequency data.
#'  \emph{ICLARM Conf. Proc}, pages 53--74, 1987.
#'   
#' @aliases powh-method powh,numeric,numeric-method
#'  
#' @export
#' @docType methods
#' @rdname powh
#' 
#' @examples
#' \dontrun{
#' data(bonLn)
#' rslt=with(subset(bonLn,year==2013), powh(len,n))
#' }
setGeneric('powh', function(len,n,...) standardGeneric('powh'))
setMethod("powh", signature(len='numeric', n="numeric"),
          function(len,n,weights=FALSE,fromMode=FALSE,linf=0){
  
  fn=function(len,n){
    #require(plyr)
    
    res=ddply(data.frame(n=n,len=len), .(len), function(x) sum(x$n))
    res=res[order(res$len),]
    
    csum =rev(cumsum(rev(res$V1)))
    clsum=rev(cumsum(rev(res$len*res$V1)))
    mn   =clsum/csum
    
    data.frame(mn=mn,diff=mn-res$len,len=res$len,n=res$V1)}
  
  linfFn=function(coef) -coef[1]/coef[2]
  zkFn  =function(coef) (-1-coef[2])/coef[2]
  
  if (fromMode){
    mode=seq(length(n))[max(n)==n]
    len=len[mode:length(len)]
    n  =  n[mode:length(n)]
  }
  
  dat=fn(len-linf,n)
  
  if (!weights)
    if (all(linf==0))
      res=lm(diff~len,data=dat)
    else
      res=lm(diff~0+len,data=dat)
  else   
    if (all(linf==0))
      res=lm(diff~len,weights=n,data=dat)
    else  
      res=lm(diff~0+len,weights=n,data=dat)
 
  if (all(linf==0))
    coef=res$coefficients
  else  
    coef=c(linf,res$coefficients)
  
  if (all(linf==0))
     params=c("linf"=linfFn(coef),"zk"=zkFn(coef))
  else
     params=c("linf"=linf,"zk"=zkFn(coef))
  
  names(params)=c("linf","zk")
  
  dat=dat[is.finite(dat$mn),]
  dat$len=dat$len+linf
  dat$mn =dat$mn +linf
  predict=predict(res,len=dat$len)

  return(list(params=params,data=data.frame(dat,hat=predict)))})


#' moment
#' 
#' @description 
#'    aa
#'      
#' @param object a vector holding a time series
#' @param ... any other arguments, i.e. x,n=rep(1,length(x)),na.rm=T
#' 
#' @return a \code{vector} with the inter-annual variation each time step
#' @export
#' @docType methods
#' @rdname utils
#' 
#' @aliases moment-method moment,numeric-method 
#' 
#' @examples
#' x=1
setGeneric('moment', function(object,...) standardGeneric('moment'))
setMethod("moment", signature(object='numeric'),
          function(object,n=rep(1,length(object)),na.rm=T) { 
  if(length(n)==1) n=rep(n,length(object)) 
  
  mn= sum(object*n,            na.rm=na.rm)/sum(n,na.rm=na.rm)
  sd=(sum(n*(object-mn)^2,     na.rm=na.rm)/sum(n,na.rm=na.rm))^.5
  sk= sum(n*((object-mn)/sd)^3,na.rm=na.rm)/sum(n,na.rm=na.rm)
  ku= sum(n*((object-mn)/sd)^4,na.rm=na.rm)/sum(n,na.rm=na.rm)-3
  
  ## weighted median
  n=unlist(c(n))
  object=unlist(c(object))
  
  cumn=cumsum(n)/sum(n)
  max.=max((1:length(n))[cumn<.50])
  min.=min((1:length(n))[cumn>=.50])
  
  me=(x[min.]*n[min.]+x[max.]*n[max.])/(n[min.]+n[max.])
  
  return(c(mn=mn,sd=sd,sk=sk,ku=ku,me=me))})

#' unbin
#' 
#' @description 
#'  For a vector with labels corresponding to intervals i.e. \code{"(0,10]"}
#'  returns a data.frame with left and right boundaries and mid point.
#'      
#' @param x a vector of with intervals as names 
#' 
#' @return a \code{data.frame} with left and right boundaries and mid points.
#' @export
#' @docType methods
#' @rdname powh
#' 
#' @examples
#' x=summary(cut(runif(100),seq(0,1,.1)))
#' unbin(x)
unbin=function(x){
  nms  =names(x)
  left =as.numeric(substr(nms,2,unlist(gregexpr(",",nms))-1))
  right=as.numeric(substr(nms,  unlist(gregexpr(",",nms))+1,nchar(nms)-1))
  mid  =(left+right)/2
  
  data.frame(left=left,right=right,mid=mid,n=x)}
