globalVariables("asym")

mf2FLPar=function(x){
  dmns=dimnames(x)[2:1]
  names(dmns)=c("params","iter")
  x=t(as.matrix(x))
  
  FLPar(array(x,dim=dim(x),dimnames=dmns),units="")}
#' lhPar
#' 
#' Uses life history theory to derive parameters for biological relationships, i.e. or growth, maturity, natural mortality
#' 
#' @description
#'     
#' Selectivity by default is set so age at peak selectivity is the same as age at 50/% mature (a50) As a minimum all `lhPar` requires is `linf` the asymptotic length of the von Bertalannfy growth equation. 
#'  
#' @param   params \code{FLPar} object with parameters for life history equations and selection pattern.
#' Need Linfinity to estimate other parameters, if any other parameters supplied in \code{code} then
#' these are not provided by the algorithm 
#' @param   t0 of von Bertalanffy. This is a default that isnt normally derived
#' from life history theory, as are the following args.
#' @param   a coefficient of length weight relationship
#' @param   b exponent of length weight relationship
#' @param   ato95 age at which 95\% of fish are mature, offset to age at which 50\% are mature
#' @param   sl selectivity-at-age parameter, standard deviation of lefthand limb of double normal
#' @param   sr stock recruitment relationship
#' @param   s steepness of stock recruitment relationship
#' @param   v virgin biomass
#' 
#' @export
#' 
#' @import methods
#' @docType methods
#' @rdname lhPar
#' @return An \code{FLPar} object with parameters
#' @examples
#' \dontrun{
#' lhPar(FLPar(linf=200))
#' }
#' 
# \deqn{ f(x) = \left\{
# \begin{array}{ll}
# 0 & x < 0 \\
# 1 & x \ge 0
# \end{array}
# \right. }{ (non-Latex version) }
lhPar=function(params,t0=-0.1,a=0.001,b=3,ato95=1,sl=2,sr=5000,s=0.9,v=1000){
 
  if("data.frame"%in%class(params)) params=mf2FLPar(params)
  
  #attach(list(t0=-0.1,a=0.001,b=3,ato95=1,sl=2,sr=5000,s=0.9,v=1000))
  
  names(dimnames(params)) <- tolower(names(dimnames(params)))
  
  if (!("t0"    %in% dimnames(params)$params)) params=rbind(params,FLPar("t0"    =t0, iter=dims(params)$iter))
  if (!("a"     %in% dimnames(params)$params)) params=rbind(params,FLPar("a"     =a,  iter=dims(params)$iter))
  if (!("b"     %in% dimnames(params)$params)) params=rbind(params,FLPar("b"     =b,  iter=dims(params)$iter))
  if (!("asym"  %in% dimnames(params)$params)) params=rbind(params,FLPar("asym"  =1,  iter=dims(params)$iter))
  if (!("bg"    %in% dimnames(params)$params)) params=rbind(params,FLPar("bg"    =b,  iter=dims(params)$iter))
  if (!("sl"    %in% dimnames(params)$params)) params=rbind(params,FLPar("sl"    =sl, iter=dims(params)$iter))
  if (!("sr"    %in% dimnames(params)$params)) params=rbind(params,FLPar("sr"    =sr, iter=dims(params)$iter))
  if (!("s"     %in% dimnames(params)$params)) params=rbind(params,FLPar("s"     =s,  iter=dims(params)$iter))
  if (!("v"     %in% dimnames(params)$params)) params=rbind(params,FLPar("v"     =v,  iter=dims(params)$iter))

  ## growth parameters
  if (!("k"     %in% dimnames(params)$params)) params=rbind(params,FLPar("k"=3.15*params["linf"]^(-0.64), iter=dims(params)$iter)) # From Gislason et al 2008, all species combined
  
  # Natural mortality parameters from Model 2, Table 1 Gislason 2010
  if (!all(c("m1","m2")%in%dimnames(params)$params))
    params=rbind(params,FLPar(m1= 0.55*(params["linf"]^1.44)%*%params["k"], iter=dims(params)$iter),
                        FLPar(m2=-1.61                                    , iter=dims(params)$iter))

  if (!("ato95" %in% dimnames(params)$params)) params=rbind(params,FLPar("ato95" =ato95, iter=dims(params)$iter))
  if (!("sl"    %in% dimnames(params)$params)) params=rbind(params,FLPar("sl"    =sl,    iter=dims(params)$iter))
  if (!("sr"    %in% dimnames(params)$params)) params=rbind(params,FLPar("sr"    =sr,    iter=dims(params)$iter))
 
  ## maturity parameters from http://www.fishbase.org/manual/FishbaseThe_MATURITY_Table.htm
  if (!("asym"    %in% dimnames(params)$params)) params=rbind(params,FLPar("asym"    =asym, iter=dims(params)$iter))

  if (!("a50" %in% dimnames(params)$params)){
    if (!("l50" %in% dimnames(params)$params)){
      l50=0.72*params["linf"]^0.93
      dimnames(l50)$params="l50"
      }else{
      l50=params["l50"]
      }
    
    a50=log(1-(l50%/%params["linf"]))%/%(-params["k"])%+%params["t0"]
    dimnames(a50)$params="a50"
    
    params=rbind(params,a50)
    }

  ## selectivity guestimate
  a1=params["a50"]
 
  dimnames(a1)$params="a1"
 
  params=rbind(params,a1)
  
  attributes(params)$units=c("cm","kg","1000s")
  
  order=c("linf","k","t0","a","b","ato95","a50","asym","bg","m1","m2","a1","sl","sr","s","v")  
  order=order[order%in%dimnames(params)[[1]]]
  order=c(order,dimnames(params)[[1]][!(dimnames(params)[[1]]%in%order)])

  return(params[order])}

setUnits=function(res, par){

    if (is.null(attributes(params)$units)) return(res)
    units=attributes(params)$units
    #browser()
    allUnits=list("params"=      "",          
               "refpts"=         "",            
               "fbar"=           "",        
               "fbar.obs"=       "",    
               "landings.obs"=   paste(units[2],units[3]),
               "discards.obs"=   paste(units[2],units[3]),
               "rec.obs"=        units[3],         
               "ssb.obs"=        paste(units[2],units[3]),
               "stock.obs"=      paste(units[2],units[3]),
               "profit.obs"=     "",     
 #              "revenue.obs"=    "",    
               "landings.sel"=   "",    
               "discards.sel"=   "", 
               "bycatch.harvest"="",        
               "stock.wt"=       units[2],     
               "landings.wt"=    units[2],     
               "discards.wt"=    units[2],      
               "bycatch.wt"=     units[2],               
               "m"=              "",             
               "mat"=            "proportion", 
               "harvest.spwn"=   "proportion",          
               "m.spwn"=         "proportion",    
               "availability"=   "proportion",           
               "price"=          "",           
               "vcost"=          "",           
               "fcost"=          "")            

  
 res@units[names(allUnits)]=allUnits
    
 return(res)}
