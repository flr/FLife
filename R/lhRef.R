#' @title Reference points based on life histories
#' 
#' @description 
#' \code{lhRef} calculates a variety of reference points i.e.  population growth rate at small 
#' population sizs (r), and at B~{MSY} (rc), ratio of virgin biomass to B~{MSY} (sk), life 
#' time reproductive output (srp0) and reproductive output at B~{MSY} (sprmsy)
 #' 
#' @param params \code{FLPar}
#' @param m function for natural mortality
#' @param sr \code{character} with stock recruitment relationship, e.g. "bevholt","ricker",...
#' @param range ages used i.e. c(min=0,max=40,minfbar=1,maxfbar=40,plusgroup=40)
#' @param what quantities to calculate "r","lopt","rc","sk","spr0","sprmsy"
#' @param msy, \code{character} with "msy", "f0.1", ...
#' 
#' @return object of class \code{FLPar} with reference points, i.e r, rc, sk, lopt,
#' 
#' @import FLRP
#' 
#' @export
#' @docType methods
#' @rdname lhRef
#' 
#' @seealso \code{\link{lhPar}}, \code{\link{lhEql}}  
#' 
#' @examples
#' \dontrun{
#' library(FLRP)
#' params=FLPar(linf=100,t0=0,k=.4)
#' params=lhPar(params)
#' lhRef(params)
#' }
lhRef<-function(params,
                 m=function(length,params) exp(0.55)*(length^-1.61)%*%(params["linf"]^1.44)%*%params["k"],
                 sr="bevholt",
                 range=c(min=0,max=40,minfbar=1,maxfbar=40,plusgroup=40),
                 what=c("r","rc","msy","lopt","sk","spr0","sprmsy"),
                 msy="msy"){
  
  eql=lhEql(params,m=m,sr=sr,range=range)

  res=NULL
  if ("r"%in%what){
    fcrash=FLQuant(refpts(eql)["crash","harvest",drop=T],
                   dimnames=list(iter=seq( dims(eql)$iter)))
    
    fcrash[is.na(fcrash)]=4
    rate =maply(dimnames(params)$iter, function(iter) 
               log(lambda(leslie(iter(eql,iter),fbar=iter(fcrash,iter))[,,1,drop=T])))

    res=cbind(res,"r"=c(rate))
    }

  if ("rc"%in%what){
    fmsy  =FLQuant(refpts(eql)[msy,"harvest",drop=T],
                 dimnames=list(iter=seq( dims(eql)$iter)))

    ratec =maply(dimnames(params)$iter, function(iter) 
      log(lambda(leslie(iter(eql,iter),fbar=iter(fmsy,iter))[,,1,drop=T])))
    
    res=cbind(res,"rc"     =c(ratec))
    }
 
  if ("msy"%in%what){
    msy. =FLQuant(refpts(eql)[msy,"yield",drop=T],
                   dimnames=list(iter=seq( dims(eql)$iter)))
    
    res=cbind(res,"msy"     =c(msy.))
  }
  
  if ("lopt"%in%what){
    lopt=loptAge(params,m=m)
    
    res=cbind(res,"lopt"   =c(lopt))
    }
  
  if ("sk"%in%what){
    sk   =refpts(eql)[msy,"ssb"]/refpts(eql)["virgin","ssb"]
    
    res=cbind(res,"sk"   =c(sk))
    }

  if ("spr0"%in%what){
   spr0=refpts(eql)["virgin","ssb"]/refpts(eql)["virgin","rec"]
  
   res=cbind(res,"spr0"    =c(spr0))}

  if ("sprmsy"%in%what){
    sprmsy=refpts(eql)[msy,"ssb"]/refpts(eql)[msy,"rec"]
    
    res=cbind(res,"sprmsy"    =c(sprmsy))}
  
  res=t(res)
  names(dimnames(res))=c("params","iter")
  FLPar(res)}


if (FALSE){

library(ggplot2)
library(FLCore)
library(FLRP)
library(FLife)
library(popbio)
library(ggplotFL)
library(RSQLite)
library(DBI)

source('~/Desktop/flr/git/FLife/R/leslie.R')
source('~/Desktop/flr/git/FLife/R/lopt.R')

drv =dbDriver("SQLite")

db  ="/home/laurie/Desktop/Dropbox/Largo/FLife-runs/data/mis3.db"
#dbDisconnect(con)
con=dbConnect(drv, dbname=db)
res=dbReadTable(con,"res")

spp=c("Cheimerius nufar","Gadus morhua","Lutjanus campechanus",       
      "Pristipomoides filamentosus","Rastrelliger kanagurta","Solea solea")

scen=expand.grid(juve =!c(TRUE,FALSE),
                 dome =!c(TRUE,FALSE),
                 steep=c(0.75,.9),
                 m    =c("Gislason","Constant"),
                 species=spp,
                 stringsAsFactors=FALSE)

load("/home/laurie/Desktop/Dropbox/Largo/FLife-runs/data/imp.RData")
imp=subset(imp,species%in%spp&Set=="Original")[,-2]
scen=merge(scen,imp,by=c("species"),all=TRUE)
scen=scen[do.call(order,scen[,c("species","juve","dome","steep","m","iter")]),]
scen=merge(scen,cbind(res[,c(1:4,6:7)],flag=1),by=c("species","juve","dome",
                                                    "m","steep","iter"),all.x=TRUE)
scen=subset(scen,is.na(flag))
table(scen$species)

tmp=d_ply(scen,.(juve,dome,steep,m,species), with, {
  
            
  Set="Original"
  l50=pmin(l50,linf*.8)
  
  p=data.frame("linf"=linf,"k"=k,"t0"=-0.1,"l50"=l50,"a"=a,"b"=b)
  p=gislasim(p)

  p["s",] =steep
  p["a1",]=ifelse(juve,.75,1)*p["a1",]
  p["sr",]=ifelse(dome,5,5000)
  
  if (m=="Gislason"){
    eql=lhEql(p)
  }else{
    eql=lhEql(p,m=function(params,length){
      length[]=params["l50"]
      0.55*(length^-1.66)%*%(params["linf"]^1.44)%*%params["k"]})}
  
  ref=refpts(eql)
  f  =ref["crash","harvest"]
  f[is.na(f)]=4
  
  lsl=leslie(eql,fbar=FLQuant(f[drop=T],dimnames=list(iter=iter)))
  rate=r(lsl)
                 
  f    =ref["f0.1","harvest"]
  ratec=r(leslie(eql,fbar=FLQuant(f[drop=T],dimnames=list(iter=iter))))
              
  if (m=="Gislason"){
     lop  =try(lopt(p))
  }else{
     lop  =lopt(p,
          m=function(length,params){
                        length[]=params["l50"]
                        0.55*(length^-1.66)%*%(params["linf"]^1.44)%*%params["k"]}
                   )}
  
  spr0=eql
  params(spr0)=params(eql)["a"]/params(eql)["a"]
  model( spr0)=geomean()$model
  spr0        =computeRefpts(spr0)
  
  res=data.frame("juve"   =juve,
                 "dome"   =dome,
                 "steep"  =steep,
                 "m"      =m,
                 "Set"    =Set,
                 "species"=species,                   
                 "iter"   =iter,
                 "r"      =c(rate),
                 "rc"     =c(ratec),
                 "lopt"   =c(lop),
                 "sk"     =c(ref["msy","ssb"]/ref["virgin","ssb"]),
                 "spr0"    =c(spr0["msy","ssb"]/spr0["virgin","ssb"]))

  dbWriteTable(con,"res",res,append=TRUE)})
  }