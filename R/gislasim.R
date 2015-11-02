mf2FLPar=function(x){
  dmns=dimnames(x)[2:1]
  names(dmns)=c("params","iter")
  x=t(as.matrix(x))
  
  FLPar(array(x,dim=dim(x),dimnames=dmns),units="")}

#' gislasim
#' 
#' Takes an \code{FLPar} object with life history and selectivity parameters
#' and generates an corresponding \code{FLBRP} object. Can uses a range of functional forms
#' 
#' @param   \code{par} an \code{FLPar} object with life history parameters
#' @param   \code{growth} function for growth
#' @param   \code{fnM} function for natutal mortality      
#' @param   \code{fnMat} function for proportion mature-at-age
#' @param   \code{fnSel} function for selectivity-at-age
#' @param   \code{sr} character for stock recruitment relationship
#' @param   \code{range} age range
#' @param   \code{spwn} proportion of year when spawning occurrs, i.e. level of natural mortality prior to spawning
#' @param   \code{fish} proportion of year when fishing happens
#' @param   \code{units} 
#' @param   ...
#' 
#' @return \code{FLBRP} 
#' 
#' @export
#' @docType methods
#' @rdname gislasim
#' 
#' @seealso \code{\link{vonB}} \code{\link{lorenzen}} \code{\link{logistic}} \code{\link{doubleNormal}}  
#' 
#' @examples
#' \dontrun{
#' par=gislasim(FLPar(linf=100))
#' }
gislasim=function(par,t0=-0.1,a=0.001,b=3,ato95=1,sl=2,sr=5000,s=0.9,v=1000){
 
  if("data.frame"%in%class(par)) par=mf2FLPar(par)
  
  #attach(list(t0=-0.1,a=0.001,b=3,ato95=1,sl=2,sr=5000,s=0.9,v=1000))
  
  names(dimnames(par)) <- tolower(names(dimnames(par)))
  
  if (!("t0"    %in% dimnames(par)$params)) par=rbind(par,FLPar("t0"    =t0, iter=dims(par)$iter))
  if (!("a"     %in% dimnames(par)$params)) par=rbind(par,FLPar("a"     =a,  iter=dims(par)$iter))
  if (!("b"     %in% dimnames(par)$params)) par=rbind(par,FLPar("b"     =b,  iter=dims(par)$iter))
  if (!("asym"  %in% dimnames(par)$params)) par=rbind(par,FLPar("asym"  =1,  iter=dims(par)$iter))
  if (!("bg"    %in% dimnames(par)$params)) par=rbind(par,FLPar("bg"    =b,  iter=dims(par)$iter))
  if (!("sl"    %in% dimnames(par)$params)) par=rbind(par,FLPar("sl"    =sl, iter=dims(par)$iter))
  if (!("sr"    %in% dimnames(par)$params)) par=rbind(par,FLPar("sr"    =sr, iter=dims(par)$iter))
  if (!("s"     %in% dimnames(par)$params)) par=rbind(par,FLPar("s"     =s,  iter=dims(par)$iter))
  if (!("v"     %in% dimnames(par)$params)) par=rbind(par,FLPar("v"     =v,  iter=dims(par)$iter))

  ## growth parameters
  if (!("k"     %in% dimnames(par)$params)) par=rbind(par,FLPar("k"=3.15*par["linf"]^(-0.64), iter=dims(par)$iter)) # From Gislason et al 2008, all species combined
  
  # Natural mortality parameters from Model 2, Table 1 Gislason 2010
  if (!all(c("m1","m2")%in%dimnames(par)$params))
    par=rbind(par,FLPar(m1= 0.55*(par["linf"]^1.44)%*%par["k"], iter=dims(par)$iter),
                  FLPar(m2=-1.61                              , iter=dims(par)$iter))

  if (!("ato95" %in% dimnames(par)$params)) par=rbind(par,FLPar("ato95" =ato95, iter=dims(par)$iter))
  if (!("sl"    %in% dimnames(par)$params)) par=rbind(par,FLPar("sl"    =sl,    iter=dims(par)$iter))
  if (!("sr"    %in% dimnames(par)$params)) par=rbind(par,FLPar("sr"    =sr,    iter=dims(par)$iter))
 
  ## maturity parameters from http://www.fishbase.org/manual/FishbaseThe_MATURITY_Table.htm
  if (!("asym"    %in% dimnames(par)$params)) par=rbind(par,FLPar("asym"    =asym, iter=dims(par)$iter))

  if (!("a50" %in% dimnames(par)$params)){
    if (!("l50" %in% dimnames(par)$params)){
      l50=0.72*par["linf"]^0.93
      dimnames(l50)$params="l50"
      }else{
      l50=par["l50"]
      }
    
    a50=log(1-(l50%/%par["linf"]))%/%(-par["k"])%+%par["t0"]
    dimnames(a50)$params="a50"
    
    par=rbind(par,a50)
    }

  ## selectivity guestimate
  a1=par["a50"]
 
  dimnames(a1)$params="a1"
 
  par=rbind(par,a1)
  
  attributes(par)$units=c("cm","kg","1000s")
  
  return(par)}

setUnits=function(res, par){

    if (is.null(attributes(par)$units)) return(res)
    units=attributes(par)$units
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