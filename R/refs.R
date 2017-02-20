#' @title refs
#' 
#' @description 
#' 
#' Returns a variety of reference points based on spawner per recruit (S/R) and stock recruitment
#' relationships.
#' 
#' $F_{MSY}$ corresponds to the level of exploitation that provides the maximum, derived from a 
#' spawner and yield curves combined with a stock recruitment relationship, while $F_{Crash}$ is 
#' the level of F that will drive the stock to extinction. 
#' 
#' Reference points can also be derived from the spawner and yield curves alone. For example
#' $F_{0.1}$ is a proxy for $F_{MSY}$ and is the fishing mortality that corresponds to a point 
#' on the yield per recruit curve where the slope is 10\% of that at the origin; $F_{Max}$ is 
#' the maximum of the yield per recruit and $SPR_O$ is the spawner per recruit at virgin biomass, 
#' and $SPR_30$ corressponds to the point on the curve where SPR is 30\% of $SPR_0$. In these 
#' cases the biomass, ssb and yield values are derived by multiplying the per recruit values by 
#' the average recruitment. 
#'
#' @param object \code{FLStock} 
#' 
#' @return \code{data.frame} 
#' 
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
  
  warn=options()$warn
  options(warn=-1)
  
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
  
  if (dims(refpts(eql))$iter>1){
    catch(object)=propagate(catch(object),dims(refpts(eql))$iter)
    dimnames(catch(object))$iter=dimnames(rec(object))$iter
    }
  
  current=as.data.frame(mcf(FLQuants(object,"harvest"=fbar,"yield"=catch,
                                            "rec"    =rec, "ssb"=ssb,
                                            "biomass"=computeStock)),drop=TRUE)[,-1]
  
  
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
            b   =subset(res,quantity=="biomass")[,-(1:2)],
            s   =subset(res,quantity=="ssb"    )[,-(1:2)],
            r   =subset(res,quantity=="rec"    )[,-(1:2)],
            f   =subset(res,quantity=="harvest")[,-(1:2)],
            y   =subset(res,quantity=="yield"  )[,-(1:2)])
  
  res=transform(merge(res,r,by="iter"),
                 rt=log(b.msy/b.current)/r)
  
  options(warn=warn)
  
  res[, !(names(res)%in%c("b.crash","s.crash","r.crash","y.crash",
                          "b.year", "s.year", "r.year", "y.year","f.year",
                          "f.virgin","f.spr.0","y.virgin","y.spr.0"))]})
