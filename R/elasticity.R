globalVariables(c("fwd","jacobian","melt"))
globalVariables(c("fwd","jacobian","melt"))
globalVariables(c("maply","jacobian","br","melt","maply","lambda","mdply","jacobian"))
    

#' elasticity
#'
#' Estimates elasticity
#'  
#' @param params parameters
#' @param sel selection pattern
#' @param fn user supplied function
#' 
#' @aliases elasticity-method
#' 
#' @return \code{FLPar} 
#' 
#' #' @export
#' @docType methods
#' @rdname elasticity
#'   
#' @seealso \code{\link{lhPar}}  
#'  
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(FLCore)
#' library(FLBRP)
#' library(FLife)
#' data(pars)
#' pms=lhSim(pars[[1]])
#' eql=lh(pms,range = c(min=0,max=8, minfbar=1,maxfbar=8,plusgroup=8))
#' lsl=leslie(eql,fbar=c(refpts(eql)["crash","harvest"]))
#' }
elasticity=function(params,sel,fn){

   elasFn=function(x,dmns,what,sel,fn) {

    pr       =exp(x)
    names(pr)=dmns
    br       =fn(FLPar(pr),sel)
      
    br       =brp(br)
    rp       =br@refpts
 
    smy=c(c(rp[what,"ssb"]),
          c(rp[what,"biomass"]),
          c(rp[what,"harvest"]),
          c(rp[what,"yield"]),
          c(rp[what,"rec"]),
          c(ssb(    br)),
          c(stock(  br)),
          c(fbar(   br)),
          c(catch(  br)),
          c(rec(    br)),
          c(ssb(    br)/rp[what,"ssb"]),
          c(stock(  br)/rp[what,"biomass"]),
          c(fbar(   br)/rp[what,"harvest"]),
          c(catch(  br)/rp[what,"yield"]),
          c(rec(    br)/rp[what,"rec"]),
          maply(c(fbar(br)), function(x,br) log(lambda(leslie(br,x))),br=br))

    return(log(smy))}

   jbn=jacobian(elasFn,log(c(par)),dmns=dimnames(par)$params,what="msy",sel=sel,fn=fn)
   
   lbl=data.frame(Year    =c(rep(1,5),rep(seq(dims(fbar(br))$year),11)),
                  Quantity=c(c("SSB","Biomass","Harvest","Yield","Recruits"),rep(rep(c("SSB","Biomass","Harvest","Yield","Recruits"),each=dims(fbar(br))$year),2),rep("r",each=dims(fbar(br))$year)),
                  Type    =c(rep("Reference Point",5),rep(c("Absolute","Relative"),each=dims(fbar(br))$year*5),rep("Pop",dims(fbar(br))$year)))
  
   res=cbind(lbl,jbn)
   res=melt(res,id.var=c("Year","Quantity","Type"))
   res$Parameter      =factor(dimnames(par)$params[res$variable])
            
   return(res[,c("Year","Quantity","Type","Parameter","value")])}

elasticity2=function(par,sel,fn){

   elasFn=function(x,dmns,sel,fn) {

    pr       =exp(x)
    names(pr)=dmns
    br       =fn(FLPar(pr),sel)
      
    br       =brp(br)
    
    smry=function(what,br){
          
        rp =br@refpts
 
        res=c(c(rp[what,"ssb"]),
              c(rp[what,"biomass"]),
              c(rp[what,"harvest"]),
              c(rp[what,"yield"]),
              c(rp[what,"rec"]),
              c(ssb(    br)),
              c(stock(  br)),
              c(fbar(   br)),
              c(catch(  br)),
              c(rec(    br)),
              c(ssb(    br)/rp[what,"ssb"]),
              c(stock(  br)/rp[what,"biomass"]),
              c(fbar(   br)/rp[what,"harvest"]),
              c(catch(  br)/rp[what,"yield"]),
              c(rec(    br)/rp[what,"rec"]),
              maply(c(fbar(br)), function(x,br) log(lambda(leslie(br,x))),br=br))
        
        log(res)}
        
      smy=mdply(dimnames(refpts(br))$refpt, smry, br=br)
          
      return(unlist(c(smy)))}

   #tmp=elasFn(log(c(par)),dmns=dimnames(par)$params,sel=sel,fn=fn)
 
   jbn=jacobian(elasFn,log(c(par)),dmns=dimnames(par)$params,sel=sel,fn=fn)
   
   lbl=data.frame(Year    =c(rep(1,5),rep(seq(dims(fbar(br))$year),11)),
                  Quantity=c(c("SSB","Biomass","Harvest","Yield","Recruits"),rep(rep(c("SSB","Biomass","Harvest","Yield","Recruits"),each=dims(fbar(br))$year),2),rep("r",each=dims(fbar(br))$year)),
                  Type    =c(rep("Reference Point",5),rep(c("Absolute","Relative"),each=dims(fbar(br))$year*5),rep("Pop",dims(fbar(br))$year)))
  
   res=cbind(lbl,jbn)
   res=melt(res,id.var=c("Year","Quantity","Type"))
   res$Parameter      =factor(dimnames(par)$params[res$variable])
            
   return(res[,c("Year","Quantity","Type","Parameter","value")])}


setMethod('sv', signature(x='FLPar', model='character'),
          function(x, model, spr0=NA){
            
            a=x["a"]
            b=x["b"]
            s=FLPar(a=1,dimnames=dimnames(a))  
            v=FLPar(b=1,dimnames=dimnames(a))  
            spr0=FLPar(spr0,dimnames=dimnames(a))  
            
            if ("spr0" %in% dimnames(x)$params)
              spr0=x["spr0"] 
            
            c=FLPar(c=1,dimnames=dimnames(a))  
            d=FLPar(d=1,dimnames=dimnames(a))  
            if (("c" %in% dimnames(x)$params))  c=x["c"]
            if (("d" %in% dimnames(x)$params))  d=x["d"]
            
            v <- v*spr2v(model, spr0, a, b, c, d)
            s <- s*srr2s(model, ssb=v*.2, a=a, b=b, c=c, d=d) / srr2s(model, ssb=v, a=a, b=b, c=c, d=d)
            
            res=rbind(s, v, spr0)
            
            if ("c" %in% dimnames(x)$params)
              res=rbind(res, c)
            
            if ("d" %in% dimnames(x)$params)
              res=rbind(res, d)
            
            res=rbind(res, spr0)
            
            return(res)})

#x=FLPar(s=0.75,v=1500,spr0=12)
#sv(ab(x,"bevholt"),"bevholt")


doIt=function(what,par,dynamic=FALSE,fbar=FLQuant(c(seq(0,.75,length.out=21),seq(.75,.0,length.out=21)[-1]),dimnames=list(year=1:41))){
  
  #require(reshape)
  
  func=function(x,dmns,par,fbar,what) {
    
    unt=units(par) 
    par[dmns] =exp(x)
    par["t0"] =-par["t0"]
    par["M2"] =-par["M2"]
    units(par)=unt  
    
    res=lh(par)
    
    fbar(res)=fbar
    res      =brp(res)
    rp       =refpts(res)
    if (dynamic) res=fwd(res)
    
    smy=c(c(rp[what,"ssb"]),
          c(rp[what,"biomass"]),
          c(rp[what,"harvest"]),
          c(rp[what,"yield"]),
          c(rp[what,"rec"]),
          c(ssb(    res)),
          c(stock(  res)),
          c(fbar(   res)),
          c(catch(  res)),
          c(rec(    res)),
          c(ssb(    res)/rp[what,"ssb"]),
          c(stock(  res)/rp[what,"biomass"]),
          c(fbar(   res)/rp[what,"harvest"]),
          c(catch(  res)/rp[what,"yield"]),
          c(rec(    res)/rp[what,"rec"]))
    
    return(log(smy))
  }
  
  jbn=jacobian(func,log(c(par)),dmns=dimnames(par)$params,par=par,fbar=fbar,what=what)
  
  res=data.frame(Year    =c(rep(1,5),rep(seq(dims(fbar)$year),10)),
                 Quantity=c(c("SSB","Biomass","Harvest","Yield","Recruits"),rep(rep(c("SSB","Biomass","Harvest","Yield","Recruits"),each=dims(fbar)$year),2)),
                 Type    =c(rep("Reference Point",5),rep(c("Absolute","Relative"),each=dims(fbar)$year*5)))
  
  res=cbind(res,jbn)
  res=melt(res,id.var=c("Year","Quantity","Type"))
  res$Parameter      =factor(dimnames(par)$params[res$variable])
  
  p.=data.frame(Parameter=c("linf",  "t0",    "M1","M2","s",  "vb", "a",     "b",     "bg",      "k",       "ato95",      "sl",         "sr",      "a50",     "asym",       "a1",         "fec"),
                Process  =c("Growth","Growth","M", "M", "SRR","SRR","Growth","Growth","Maturity","Growth",  "Maturity","Selectivity","Selectivity","Maturity","Maturity","Selectivity","Maturity"))
  
  res=merge(res,p.)
  
  return(res[,c("Year","Quantity","Type","Parameter","Process","value")])}

cvIt=function(what,par,dynamic=FALSE,fbar=FLQuant(c(seq(0,.75,length.out=21),seq(.75,.0,length.out=21)[-1]),dimnames=list(year=1:41))){
  
  func=function(x,par,fbar,dynamic=TRUE,what) {
    
    par[] =x
    res=lh(par)
    
    fbar(res)=fbar
    res      =brp(res)
    rp       =refpts(res)
    if (dynamic) res=fwd(res)
    
    smy=c(c(rp[what,"ssb"]),
          c(rp[what,"biomass"]),
          c(rp[what,"harvest"]),
          c(rp[what,"yield"]),
          c(rp[what,"rec"]),
          c(ssb(    res)),
          c(stock(  res)),
          c(fbar(   res)),
          c(catch(  res)),
          c(rec(    res)),
          c(ssb(    res)/rp[what,"ssb"]),
          c(stock(  res)/rp[what,"biomass"]),
          c(fbar(   res)/rp[what,"harvest"]),
          c(catch(  res)/rp[what,"yield"]),
          c(rec(    res)/rp[what,"rec"]))
    
    return(smy)}
  
  jbn=jacobian(func,c(par),par=par,fbar=fbar,what=what)
  
  res=data.frame(Year    =c(rep(1,5),rep(seq(dims(fbar)$year),10)),
                 Quantity=c(c("SSB","Biomass","Harvest","Yield","Recruits"),rep(rep(c("SSB","Biomass","Harvest","Yield","Recruits"),each=dims(fbar)$year),2)),
                 Type    =c(rep("Reference Point",5),rep(c("Absolute","Relative"),each=dims(fbar)$year*5)))
  
  res=cbind(res,jbn)
  res=melt(res,id.var=c("Year","Quantity","Type"))
  res$Parameter      =factor(dimnames(par)$params[res$variable])
  
  p.=data.frame(Parameter=c("linf",  "t0",    "M1","M2","s",  "v", "a",     "b",     "bg",      "k",       "ato95",      "sl",         "sr",      "a50",     "asym",       "a1"      ),
                Process  =c("Growth","Growth","M", "M", "SRR","SRR","Growth","Growth","Maturity","Growth",  "Maturity","Selectivity","Selectivity","Maturity","Maturity","Selectivity"))
  
  res=merge(res,p.)
  
  return(res[,c("Year","Quantity","Type","Parameter","Process","value")])}
