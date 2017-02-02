#' @title refs
#' 
#' @description 
#'
#' @param object \code{FLStock} 
#' 
#' @return \code{data.frame} 
#' ' 
#' @export
#' @rdname refs
#' 
#' @examples
#' \dontrun{
#' 
#'  library(FLCore)
#'  library(FLBRP)
#'  library(ggplotFL)
#'  library(FLife)
#'  library(popbio)
#'  
#'  data(ple4)
#'  
#'  refs(ple4)
#' } 

setGeneric('refs', function(object,...) standardGeneric('refs'))
setMethod("refs", signature(object="FLStock"),
    function(object,s=0.9){
  
  ## msy
  eql =FLBRP(object)
  spr0=mean(spr0(FLBRP(object)))
  srr =as.FLSR(object,model="bevholtSV")
  upper(srr)[1:2]=1e12
  srr=fmle(srr,
           fixed=list(s=s,spr0=spr0),
           control=list(silent=TRUE),
           method="Brent")
  model(eql) =bevholt()$model
  params(eql)=ab(params(srr),"bevholt")[c("a","b")]
  eql        =brp(eql)
  
  ## F0.1
  eql1=brp(FLBRP(object))
  params(eql1)=propagate(params(eql1),dim(refpts(eql))[3])
  params(eql1)[]=unlist(c(ddply(rod(rec.obs(eql1)),.(iter), function(x)
    mean(subset(x,regime==max(as.numeric(regime)))$data))["V1"]))
  refpts(eql1)=refpts(eql1)[c("f0.1","spr.30","virgin"),]
  dimnames(refpts(eql1))[[1]][3]="spr.0"
  eql1=brp(eql1)
  
  res =model.frame(refpts(eql )[c("msy","crash","virgin"), c("harvest","yield","rec","ssb","biomass")])
  res1=model.frame(refpts(eql1)[c("f0.1","spr.30","spr.0"),c("harvest","yield","rec","ssb","biomass")])
  res=cbind(res[,4:5],res[1:3],res1[,1:3])
  
  current=as.data.frame(mcf(FLQuants(object,"harvest"=fbar,"yield"=catch,
                                            "rec"    =rec, "ssb"=ssb,
                                            "biomass"=computeStock)),drop=TRUE)
  current=subset(current,year==max(year)&!is.na(data))[,-1]
  names(current)[names(current)=="data"]="current"
  names(current)[names(current)=="qname"]="quantity"
  if (!("iter"%in%names(current))) 
    current=cbind(iter=1,current)
  
  res=merge(res,current,by=c("iter","quantity"))
  
  r =log(aaply(leslie(eql, f=c(refpts(eql)["crash","harvest"])),3,
                function(x) lambda(x[drop=T])))
  rc=log(aaply(leslie(eql, f=c(refpts(eql)["msy","harvest"])),3,
                 function(x)   lambda(x[drop=T])))
  r =data.frame(r=r, rc=rc, iter=seq(dim(stock.n(object))[6]))
  
  res=cbind(iter=subset(res,quantity=="ssb")[ ,  1 ],
            b=subset(res,quantity=="biomass")[,-(1:2)],
            s=subset(res,quantity=="ssb"    )[,-(1:2)],
            r=subset(res,quantity=="rec"    )[,-(1:2)],
            f=subset(res,quantity=="harvest")[,-(1:2)],
            y=subset(res,quantity=="yield"  )[,-(1:2)])
  
  res=transform(merge(res,r,by="iter"),
                 rt=log(b.msy/b.current)/r)
  
  res[, !(names(res)%in%c("b.crash","s.crash","r.crash","y.crash",
                          "b.year", "s.year", "r.year", "y.year","f.year",
                          "f.virgin","f.spr.0","y.virgin","y.spr.0"))]})
