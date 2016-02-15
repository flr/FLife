#' popPar
#'
#' Von Bertalanffy growth equation
#' 
#' @param params \code{FLPar}
#' @param m function for natural mortality
#' 
#' @return Depends on the value of \code{data} 
#' 
#' @export
#' @docType methods
#' @rdname popPar
#' 
#' @seealso \code{\link{gascuel}}  
#' 
#' @examples
#' \dontrun{
#' params=FLPar(linf=100,t0=0,k=.4)
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' len=vonB(age,params)
#' age=vonB(params,length=len)
#' }
popPar<-function(params,
                  m=function(length,params) 
                         params["m1"]%*%(exp(log(length)%*%params["m2"])),
                 range=c(min=0,max=40,minfbar=1,maxfbar=40,plusgroup=40)){
  
  eql=lhEql(params,m=m,range=range)
  
  fcrash=FLQuant(refpts(eql)["crash","harvest",drop=T],
                 dimnames=list(iter=seq( dims(eql)$iter)))
  fcrash[is.na(fcrash)]=4
  f0.1  =FLQuant(refpts(eql)["f0.1","harvest",drop=T],
                 dimnames=list(iter=seq( dims(eql)$iter)))

  rate =r(leslie(eql,fbar=fcrash))
  ratec=r(leslie(eql,fbar=f0.1))
  lop  =loptAge(params,m=m)
  sk   =refpts(eql)["msy","ssb"]/refpts(eql)["virgin","ssb"]
  
  rfs=refpts(eql)[c("f0.1","virgin")]
  dimnames(rfs)[[1]][1]="ref"
  rfs["ref",-1]=NA
  lro=eql
  params(lro)=params(eql)["a"]/params(eql)["a"]
  model( lro)=geomean()$model
  refpts(lro)=rfs
  lro        =computeRefpts(lro)
  lro        =c(lro["ref","ssb"]/lro["virgin","ssb"])
            
  data.frame("r"      =c(rate),
             "rc"     =c(ratec),
             "lopt"   =c(lop),
             "sk"     =c(sk),
             "lro"    =c(lro))
}


if (FALSE){

library(ggplot2)
library(FLCore)
library(FLBRP)
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
    eql=lh(p)
  }else{
    eql=lh(p,m=function(params,length){
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
  
  lro=eql
  params(lro)=params(eql)["a"]/params(eql)["a"]
  model( lro)=geomean()$model
  lro        =computeRefpts(lro)
  
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
                 "lro"    =c(lro["msy","ssb"]/lro["virgin","ssb"]))

  dbWriteTable(con,"res",res,append=TRUE)})
  }