ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}

### Age based ref pts
## virgin, Bmsy, msy, Fmsy, Hmsy, Fcrash, shape
### Biomass based parameters
## r, rc, k, p 
benchFn<-function(object,window=5){
  
  options(warn=-1)
  goodIts=seq(dim(m(object))[2]-(window%/%2))
  if (any((1:(window%/%2))>0))
    goodIts=goodIts[-(1:(window%/%2))]
  
  ## set up FLBRP object with iters for moving average across years
  n  =dims(object)$year
  sr =fmle(as.FLSR(object,model="bevholt"),control=list(silent=TRUE))
  bry=FLBRP(object,sr=sr)
  
  slt=getSlots("FLBRP")
  for(s in names(slt[slt%in%"FLQuant"])[c(9:20)])
    slot(bry,s)=propagate(slot(bry,s),n)
  
  for(s in names(slt[slt%in%"FLQuant"])[c(12:14,16:19)])
    slot(bry,s)[]=slot(object,s)
  
  params(bry)     =propagate(params(bry),n)
  if (window>1)
    params(bry)["a"]=params(bry)["a"]*ma(c(rec.obs(bry)),n=window)/mean(rec.obs(bry))
  
  landings.sel(bry)[]=harvest(object)%*%landings.n(object)%/%catch.n(object)
  discards.sel(bry)[]=harvest(object)%*%discards.n(object)%/%catch.n(object)
  
  ref=computeRefpts(bry)[c("virgin","msy","crash"),c("ssb","biomass","yield","harvest")]
  
  res=data.frame(
    virgin=c(ref["virgin","ssb"]),
    k     =c(ref["virgin","biomass"]),
    bmsy  =c(ref[   "msy","ssb"]),
    msy   =c(ref[   "msy","yield"]),
    fmsy  =c(ref[   "msy","harvest"]),
    hmsy  =c(ref[   "msy","yield"]/ref["msy","biomass"]),
    fcrash=c(ref[ "crash","harvest"]),
    bmsyk =c(ref["msy","ssb"]/ref["virgin","ssb"]))[goodIts,]
  
  ek    =apply(stock.n(bry)[,1]%*%stock.wt(bry)[,1]%*%catch.sel(bry)%/%
                 apply(catch.sel(bry),c(2,6),max),6,sum)[,,,,,goodIts]
  
  r=mdply(goodIts,function(i)
    lambda(leslie(iter(bry,i),c(refpts(bry)["crash","harvest",i]),exploitable=!TRUE)[drop=T])-1)
  rc=mdply(goodIts,function(i)
    lambda(leslie(iter(bry,i),c(refpts(bry)["msy","harvest",i]),exploitable=TRUE)[drop=T])-1)
  
  res=cbind(res,r=r[,"V1"],rc=rc[,"V1"],ek=c(ek))
  
  fn=function(p,msy,r,k) {
    res=(msy/(r*k)-(1/(1+p))^(1/p+1))^2
    res}
  
  p=mdply(goodIts,function(i) 
    optimise(fn, c(1e-12,1e12),    
             msy=c(res[i,"msy"]),
             r  =c(res[i,"r"]),
             k  =c(res[i,"ek"]))$minimum)[,"V1"]
  
  res=cbind(year=as.numeric(ac(dimnames(catch.obs(bry))$year[goodIts])),res,p=p)
  
  options(warn=1)
  
  res}

bench<-function(object,window=5){
  
  n  =dims(object)$year
  sr =fmle(as.FLSR(object,model="bevholt"),control=list(silent=TRUE))
  bry=FLBRP(object,sr=sr)
  
  slt=getSlots("FLBRP")
  for(s in names(slt[slt%in%"FLQuant"])[c(9:20)])
    slot(bry,s)=propagate(slot(bry,s),n)
  
  for(s in names(slt[slt%in%"FLQuant"])[c(12:14,16:19)])
    slot(bry,s)[]=slot(object,s)
  
  params(bry)     =propagate(params(bry),n)
  if (window>1)
    params(bry)["a"]=params(bry)["a"]*ma(c(rec.obs(bry)),n=window)/mean(rec.obs(bry))
  
  landings.sel(bry)[]=harvest(object)%*%landings.n(object)%/%catch.n(object)
  discards.sel(bry)[]=harvest(object)%*%discards.n(object)%/%catch.n(object)
  
  tmp=model.frame(computeRefpts(bry)["msy",   c(1:2,4:5)])
  tm2=model.frame(computeRefpts(bry)["msy","ssb"]/computeRefpts(bry)["virgin","ssb"])
  tm2[,"quantity"]="lro"
  tmp=rbind.fill(tmp,tm2,
                 transform(as.data.frame(fapexAge(object),drop=T),
                           msy=data,iter=year-1929,quantity="Fapex Age")[,c(3:5)])
  tmp=rbind.fill(tmp,
                 transform(as.data.frame(computeRefpts(bry)["msy","biomass"]/
                                           computeRefpts(bry)["virgin","biomass"]),
                           msy=data,quantity="shape")[,c(2:3,5)])
  tmp=rbind.fill(tmp,
                 transform(as.data.frame(computeRefpts(bry)["msy","yield"]/
                                           computeRefpts(bry)["msy","biomass"]),
                           msy=data,quantity="hrate")[,c(2:3,5)])
  tmp=subset(tmp,iter>window%/%2&iter<n-window%/%2)
  
  bdt=maply(seq(window%/%2+1,n-window%/%2,1),function(i) 
    FLBRP2biodyn(iter(bry,i))@params[c("r","k","p")])
  
  rc=mdply(data.frame(iter=seq(window%/%2+1,n-window%/%2,1)),function(iter){ 
    data.frame(quantity="rc",msy=lambda(leslie(iter(bry,iter),
                                               f=c(refpts(iter(bry,iter))["msy","harvest"]))[drop=TRUE])-1)
  })
  
  bdt=melt(bdt)
  names(bdt)=c("iter","quantity","msy")
  bdt=rbind(tmp,bdt,rc)
  names(bdt)[1]="data"
  
  
  bdt=transform(bdt,quantity=factor(quantity,
                                    levels=c("yield","harvest","hrate","Fapex Age",
                                             "ssb","biomass","lro",
                                             "r","k","p","shape","rc"),
                                    labels=c("MSY","Fmsy","Hmsy","FApex Age",
                                             "SSBmsy","Bmsy","Lmsy",
                                             "r","k","p","shape","rc")))
  
  
  transform(bdt,year=as.numeric(ac(iter))+dims(catch(object))$minyear-1)[,-3]}

