
## Adds the ICES PA & MSY reference points to refpts and fits a SRRR
addRefs<-function(x,refs){
  x=FLPar(NA,dimnames=list(refpt=c("virgin","msy","crash",dimnames(refs)$params),
                           quant=c("harvest","yield","rec","ssb","biomass","revenue","cost","profit"),iter=1))
  
  x[maply(dimnames(x)$refpt,function(x) gregexpr("B",x)[[1]][1]==1),"ssb"]=
    refs[maply(dimnames(refs)[[1]],function(x) gregexpr("B",x)[[1]][1]>0),]
  x[maply(dimnames(x)$refpt,function(x) gregexpr("F",x)[[1]][1]==1),"harvest"]=
    refs[maply(dimnames(refs)[[1]],function(x) gregexpr("F",x)[[1]][1]>0),]
  x} 

icesRefpts<-function(x,refs=NULL,model="bevholtSV",steepness=0.7,nyears=3) {
  eq=FLBRP(x,nyears=nyears)
  
  if (gregexpr("SV",model)[[1]][1]>0){
    sr=as.FLSR(x,model=model)
    sr=fmle(sr,
            fixed=list(s=steepness,spr0=spr0(eq)),
            control=list(silent=TRUE),
            method="Brent",
            lower=c(0.001),upper=max(ssb(sr))*10)
    params(eq)=ab(params(sr),substr(model,1,gregexpr("SV",model)[[1]][1]-1))[-dim(params(sr))[1]]
    model( eq)=do.call(substr(model,1,gregexpr("SV",model)[[1]][1]-1), list())$model
    refpts(eq)=computeRefpts(eq)
  }else{
    sr=fmle(as.FLSR(x,model=model),control=list(silent=TRUE))
    params(eq)=params(sr)
    model( eq)=do.call(model, list())$model}
  
  fbar(eq)[]=seq(0,100,1)/100*c(computeRefpts(eq)["crash","harvest"])
  
  if (!is.null(refs))
    refpts(eq)=addRefs(refpts(eq),refs)
  
  eq=brp(eq)
  eq}

## Calculates the surplus production and expected yield etc for the estinates of SSB and biomass
surplusPrd<-function(x){
  nms=dimnames(refpts(x))
  nms$refpt=paste("ssb",dimnames(ssb.obs(x))$year,sep="")
  
  rfs=FLPar(array(NA,laply(nms,length),dimnames=nms))
  rfs[,"ssb",]=ssb.obs(x)
  refpts(x)=rfs
  rtn=computeRefpts(x)
  
  rtn=alply(rtn,2,FLQuant,dimnames=dimnames(ssb.obs(x)))
  names(rtn)=as.character(unlist(attributes(rtn)$split_labels))
  
  discards.obs(x)[is.na(discards.obs(x))]=0
  
  rtn$spSSB=ssb.obs(x)[,-1]-ssb.obs(x)[,-dim(ssb.obs(x))[2]]+catch.obs(x)[,-dim(ssb.obs(x))[2]]
  rtn$spBiomass=biomass.obs(x)[,-1]-biomass.obs(x)[,-dim(biomass.obs(x))[2]]+catch.obs(x)[,-dim(biomass.obs(x))[2]]
  
  rtn=mcf(as(rtn,"FLQuants"))
  
  rtn[["peSSB"]]=rtn$spSSB-rtn$yield
  rtn[["peBiomass"]]=rtn$spBiomass-rtn$yield
  
  rtn}

## Set M, wt, sel & mat to vary by iter based on annual values to look at non-stationarity
nonStn<-function(x,sr=NULL,nyears=dim(x)[2],slots=c("m","mat","stock.wt","catch.wt","catch.sel")){
  
  if (is.null(sr)) eq=FLBRP(x,nyears=nyears)
  else             eq=FLBRP(x,sr-sr,nyears=nyears)
  
  eq=propagate(eq,dim(ssb(x))[2])
  
  q2p<-function(x){
    tmp=as.data.frame(x,drop=T)
    names(tmp)[names(tmp)=="year"]="iter"
    as.FLQuant(tmp)}
  
  if ("m"%in%slots)   m(eq)=q2p(m(x))
  if ("mat"%in%slots) mat(eq)=q2p(mat(x))
  
  if ("stock.wt"%in%slots) stock.wt(   eq)=q2p(stock.wt(   x))
  if ("catch.wt"%in%slots|"landinhgs.wt"%in%slots) landings.wt(eq)=q2p(landings.wt(x))
  if ("catch.wt"%in%slots|"discards.wt" %in%slots) discards.wt(eq)=q2p(discards.wt(x))
  
  if ("catch.sel"%in%slots|"landinhgs.sel"%in%slots) landings.sel(eq)=q2p(catch.sel(x)%*%landings.n(x)%/%catch.n(x))
  if ("catch.sel"%in%slots|"discards.sel"%in%slots)  discards.sel(eq)=q2p(catch.sel(x)%*%discards.n(x)%/%catch.n(x))
  
  nms=dimnames(refpts(eq))
  nms[[1]]=c(nms[[1]],"spr.100")
  refpts(eq)=FLPar(array(NA,lapply(nms,length),nms))
  
  rtn=computeRefpts(eq)
  
  transform(as.data.frame(rtn),year=as.numeric(dimnames(x)$year[iter]))[,-3]}