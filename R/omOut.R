omStock<-function(object){
  ageSwt<-function(object) apply(stock.n(object)%*%ages(stock.n(object)),2:6,sum)%/%
                           apply(stock.n(object),2:6,sum)
  ageCwt<-function(object) apply(catch.n(object)%*%ages(catch.n(object)),2:6,sum)%/%
                           apply(catch.n(object),2:6,sum) 
  hvt   <-function(object) catch(object)/stock(object)
  
  recs  <-function(object) {res=rec(object)
                            dimnames(res)[[1]]="all"
                            res}

  res=FLQuants(object,"ssb"=ssb,"stock"=stock,"rec"=recs,"catch"=catch,"fbar"=fbar,
                      "ageSwt"=ageSwt,"ageCwt"=ageCwt,"hrate"=hvt)
  
  model.frame(mcf(res),drop=TRUE)}

omRefs<-function(object){
  
  refs=rbind(as.data.frame(object["crash",c("harvest")]),
             as.data.frame(object["virgin",c("rec","ssb")]),
             as.data.frame(object["msy",c("yield","ssb","harvest")]))
  refs=cast(refs,iter~refpt+quant,value="data")

  refs}

omSmry<-function(x,y="missing"){
  
  res=omStock(x)
  
  if ("FLBRP" %in% is(y))
    res=merge(res,omRefs(refpts(y)))
  else if ("FLPar" %in% is(y))
    res=merge(res,omRefs(y))
  
  return(res)}
