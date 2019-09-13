extendPG<-function(object,pg=40){
    
  sr=fmle(as.FLSR(object,model="geomean"),control=list(silent=TRUE))
  bpg=setPlusGroup(object,pg)
  bpg=fwd(bpg,catch=catch(bpg)[,-1],sr=sr,sr.residuals=exp(residuals(sr)))
  
  par=FLPar(c(linf=318.9,k=0.093,t0=-0.97-(6-6)/12,a=0.0000350801,b=2.878450899))
  
  catch.wt(bpg)[11:pg,]=len2wt(vonB(ages(catch.wt(bpg)),par),par[c("a","b")])[11:pg,]
  stock.wt(bpg)[11:pg,]=len2wt(vonB(ages(stock.wt(bpg)),par),par[c("a","b")])[11:pg,]
  landings.wt(bpg)[11:pg,]=len2wt(vonB(ages(catch.wt(bpg)),par),par[c("a","b")])[11:pg,]
  discards.wt(bpg)[11:pg,]=len2wt(vonB(ages(catch.wt(bpg)),par),par[c("a","b")])[11:pg,]
  
  stock(bpg)=computeStock(bpg)
  sop=bpg
  catch(bpg)=computeCatch(bpg)
  
  bpg}

