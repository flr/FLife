omStock<-function(object){
  sage<-function(object) apply(stock.n(object)%*%ages(stock.n(object)),2:6,sum)%/%
                           apply(stock.n(object),2:6,sum)
  cage<-function(object) apply(catch.n(object)%*%ages(catch.n(object)),2:6,sum)%/%
                           apply(catch.n(object),2:6,sum) 
  swt<-function(object) apply(stock.n(object)%*%stock.wt(object),2:6,sum)%/%
    apply(stock.n(object),2:6,sum)
  cwt<-function(object) apply(catch.n(object)%*%catch.wt(object),2:6,sum)%/%
    apply(catch.n(object),2:6,sum) 
  hvt   <-function(object) catch(object)/stock(object)
  
  recs  <-function(object) {res=rec(object)
                            dimnames(res)[[1]]="all"
                            res}
  
  catchJuv<-function(object) apply(catch.n(object)%*%(1-mat(object))%*%catch.wt(object),2:6,sum)

  res=FLQuants(object,"ssb"=ssb,"stock"=stock,"rec"=recs,"catch"=catch,"fbar"=fbar,
                      "swt"=swt,"cwt"=cwt,"sage"=sage,"cage"=cage,"catchJuv"=catchJuv)
  
  model.frame(mcf(res),drop=TRUE)}

omRefs<-function(object){
  
  refs=rbind(as.data.frame(object["crash",c("harvest")]),
             as.data.frame(object["virgin",c("rec","ssb")]),
             as.data.frame(object["msy",c("yield","ssb","harvest")]))
  refs=cast(refs,iter~refpt+quant,value="data")

  refs}


omSmry<-function(x,y="missing",z="missing"){
  
  res=omStock(x)
  
  if ("FLBRP" %in% is(y))
    res=merge(res,omRefs(refpts(y)))
  else if ("FLPar" %in% is(y))  
    res=merge(res,omRefs(y))
  
  if ("FLPar" %in% is(z))
    if (all(c("a","b") %in% dimnames(z)$params))
      res=merge(res,lenFn(x,z))
  return(res)}

lenFn<-function(x,y){
  sln<-function(object) apply(stock.n(object)%*%exp(log(stock.wt(object)%/%y["a"])%/%y["b"]),2:6,sum)%/%
    apply(stock.n(object),2:6,sum)
  cln<-function(object) apply(catch.n(object)%*%exp(log(catch.wt(object)%/%y["a"])%/%y["b"]),2:6,sum)%/%
    apply(catch.n(object),2:6,sum) 
  
  model.frame(FLQuants(x,"sln"=sln,"cln"=cln),drop=TRUE)}
