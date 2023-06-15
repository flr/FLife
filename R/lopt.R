utils::globalVariables(c("growth","s","FLife"))

#' @title Lopt2
#'
#' @description Lopt, the length at which a cohort achives its maximum biomass, can be used as a
#' reference point to identify growth over- or underfishing. Since taking fish below or above 
#' this size results in potential loss of yield. The total biomass of a cohort changes through
#' time as a result of gains due to an increase in mean size-at-age and losses due to natural 
#' mortality. Lopt can therefore be estimated from the natural mortality and weight-at-age vectors.
#' 
#' @param params an \code{FLPar} object with parameter values for the natural mortality and growth 
#' functions, and the exponent \code{b} of the length/weight relationship.
#' @param mFn natural mortality function, by default Gislason
#' @param growth length or weight-at-age function, by default von Bertalanffy
#' @param ... any other arguments
#'  
#' @aliases lopt lopt-method lopt,FLPar-method
#' 
#' @return \code{FLPar} with $L_{opt}$ the length at which a cohort achives its maximum biomass
#' 
#' @details Lopt is a function of growth and natural mortality-at-age and there are several 
#' approximations such as \eqn{2/3 L_{\infty}} and \eqn{L_{\infty}\frac{3}{3+k/m}}. If the life
#' history parameters and relationships are known then $L_{opt}$ can be found by finding the 
#' time (t) and hence length at which the maximum biomass is achieved i.e.
#'  \eqn{L(T)^a e^{\int_0^T m(t)}}
#' where \eqn{m(t)} can be found from the relationship of mortality at length using the relationship 
#' of Gislason, assuming the von Bertalanffy growth curve. 
#' 
#' @export
#' @docType methods
#' @rdname lopt2
#' 
#' @seealso \code{\link{gislason}}, \code{\link{vonB}}, \code{\link{lhRef}}, \code{\link{lhPar}}, \code{\link{lhEql}},  
#' 
#' @examples
#' \dontrun{
#' params=lhPar(FLPar(linf=100,k=0.1,t0=-0.1,b=3))
#' lopt(params)
#' }
setMethod("lopt2", signature(params="FLPar"),
       function(params,
                mFn=function(length,params) exp(0.55)*(length^-1.61)%*%(params["linf"]^1.44)%*%params["k"],
                growth=FLife::vonB,
                ...){   
            dmns=dimnames(params)
            dmns$params="lopt"
            dm  =dim(params)

            res=aaply(params,seq(length(dm))[-1],function(x){
                   x.=FLPar(x)

                   rtn=try(optimise(loptFn,c(.01,c(x["linf"])*.99),params=x.,maximum=TRUE,mFn=mFn)$maximum)
   
                   if ("character" %in% mode(rtn)) rtn=NA
                   rtn})
            
            FLPar(array(res,dim=c(1, dm[-1]),dimnames=dmns))})

loptFn=function(x,params,
                mFn=function(length,params) exp(0.55)*(length^-1.61)%*%(params["linf"]^1.44)%*%params["k"],
                age=0:200,
                growth=vonB){
  
  length=growth(FLQuant(age,dimnames=list(age=age)),params)
  m.    =FLQuant(mFn(length,params), dimnames=list(age=age))
  mCum  =FLQuant(aaply(m.,2,cumsum),dimnames=dimnames(m.))

  a =qmax(growth(params=params,length=FLQuant(x)),min(age))
  a =qmin(a,max(age))
 
  aMin=floor(a)
  if (is.na(aMin)) return (as.numeric(NA))
  
  aMax=ceiling(a)
  if (is.na(aMax)) return (as.numeric(NA))

  m1=mFn(vonB(aMin,params),params)
  m2=mFn(vonB(aMax,params),params)

  slope    =(m1-m2)/(aMin-aMax)
  intercept=m1-slope*aMin
  m.=(intercept+slope*a)*(a-aMin)
    
  n =exp(-mCum[ac(aMin)]-m.)
  
  c(n%*%len2wt(x,params))
  }

if (FALSE){
 library(FLife)
  
 par=readRDS("/home/laurence-kell/Desktop/projects/mydas/data/simon/pol_par.rds")
 par=par[c("linf","k","t0","a","b")]
 par["k"]=seq(0.1,0.5,length.out=12)
 par["t0"]=-0.1
 
 loFLife=lopt2(par)

 par=lhPar(par)
 eql=lhEql(par,fbar=FLQuant(0)) # ,m=function(object,param) {res=FLQuant(m(object)); res%=%0.2;res} )
 #m(eql)[]=0.2
 res=model.frame(FLQuants(eql,biomass=function(x) catch.wt(x)%*%stock.n(x),
                              wt     =function(x) catch.wt(x)),drop=TRUE)
 
 res=transform(res,iter=factor(iter,levels=as.numeric(as.character(unique(res$iter)))))
 
 res2=ddply(res,.(iter), with, 
            data.frame(biomass=max(biomass),
                       wt     =wt[biomass==max(biomass)]))
 
 wo=as(transmute(res2,lopt=wt,iter=as.numeric(as.character(iter))),"FLPar")
 par=rbind(par,lopt=loFLife,FLPar(lopt2=wt2len(c(wo),par)))
 
 ggpairs(model.frame(par[c("k","lopt","lopt2")])[,-3])
 
 
 dat=as.data.frame(stock.n(eql)%*%stock.wt(eql))
 dat$iter=as.numeric(as.character(dat$iter))
 dat$k=c(par["k",dat$iter])
 
 ggplot(dat)+
   geom_line(aes(age,data,group=k))+
   facet_wrap(~k,scale="free")


}
  
#params=par[,10]
#params["t0"]=-.1
#loptFn(x,params)

loptFn2<-function(params) params["linf"]*3/(3+exp(params["m2"])/params["k"])

#' @title Age at maximum biomass
#'
#' @description Finds length at maximum biomass
#' 
#' @param params FLPar
#' @param m A function, i.e. gislason
#' @param growth A function, i.e. vonB
#' 
#' @param ... any other arguments
#' 
#' @aliases loptAge loptAge-method loptAge,FLPar-method
#' 
#' @return \code{FLPar} with length at maximum biomass of a cohort 
#' 
#' @details There are several ways to calculate \eqn{L_{opt}}, i.e.
#' i) \eqn{{2/3}^{rds}  L_{\infty}}
#' ii) \eqn{L_{\infty}\frac{3}{3+k/m}}
#' iii) by maximising the biomass of
#' iv) from an FLBRP object by fishing at F=0 and finding age where biomass is a maximum
#' 
#' @export
#' @docType methods
#' @rdname loptAge
#' 
#' @seealso \code{\link{loptAge}}, \code{\link{lhRef}}, \code{\link{lhPar}}, \code{\link{lhEql}},  
#' 
#' @examples
#' \dontrun{
#' params=lhPar(FLPar(linf=100))
#' loptAge(params)
#' }
setMethod("loptAge", signature(params="FLPar"),
    function(params,
                   m     =function(length,params) params["m1"]%*%(exp(log(length)%*%params["m2"])),
                   growth=vonB,
                   ...){   

      loptFn=function(x,params,m){
        
        age   =0:ceiling(x)
        dmns  =list(age=age)
        length=vonB(age=FLQuant(pmin(age+0.5,x),dimnames=dmns),params=params)
        m.    =FLQuant(mFn(length,params),    dimnames=dmns)
        mCum  =FLQuant(aaply(m.,6,sum))
        n     =exp(-mCum)
        c(n*FLife::len2wt(length[ac(ceiling(x))],params))}
           
      dmns=dimnames(params)
      dmns$params="lopt"
      dm  =dim(params)
        
      res=aaply(params,seq(length(dm))[-1],function(x){
            x.=FLPar(x)
            rtn=try(optimise(loptFn,c(0,40),params=x.,maximum=TRUE,mFn=mFn)$maximum)
            if ("character" %in% mode(rtn)) rtn=NA
            rtn})
     
      vonB(as(FLPar(array(res,dim=c(1, dm[-1]),dimnames=dmns)),"FLQuant"),params)
      })


setMethod("genTime", signature(params="FLPar"),
          function(params,
                   mFn   =function(length,params) params["m1"]%*%(exp(log(length)%*%params["m2"])),
                   growth=vonB,
                   ...){   
            
            loptFn=function(x,params,m){
              
              age   =0:ceiling(x)
              dmns  =list(age=age)
              length=vonB(age=FLQuant(pmin(age+0.5,x),dimnames=dmns),params=params)
              m.    =FLQuant(mFn(length,params),    dimnames=dmns)
              mCum  =FLQuant(aaply(m.,6,sum))
              n     =exp(-mCum)
              c(n*FLife::len2wt(length[ac(ceiling(x))],params))}
            
            dmns=dimnames(params)
            dmns$params="lopt"
            dm  =dim(params)
            
            res=aaply(params,seq(length(dm))[-1],function(x){
              x.=FLPar(x)
              rtn=try(optimise(loptFn,c(0,40),params=x.,maximum=TRUE,mFn=mFn)$maximum)
              if ("character" %in% mode(rtn)) rtn=NA
              rtn})
            

            res}
          
          )

if (FALSE){
library(FLife)
library(FLBRP)
#source('~/Desktop/flr/FLife/R/vonB.R')
#source('~/Desktop/flr/FLife/R/dnormal.R')
#source('~/Desktop/flr/FLife/R/len2wt.R')
#source('~/Desktop/flr/FLife/R/lopt.R')

#Med swordfish
par=FLPar(c(linf=238.59000,k=0.18500,t0=-1.40400,a=0.00000176,b=3.37800,l50=142.00000,
            m1=0.62774,m2=-0.30900))  

pars=FLPars("Lorenzen"=lhPar(par))
eql=FLBRPs(list("Lorenzen"=
                  lhEql(pars[["Lorenzen"]],range=c(min=0,max=120,minfbar=0,maxfbar=120,plusgroup=120))))
pars[["Lorenzen"]]["m1"]=pars[["Lorenzen"]]["m1"]*0.2/m(eql[["Lorenzen"]])["3"]
pars[["0.2"]]=pars[["Lorenzen"]]
pars[["0.2"]][c("m1","m2")]=c(0.2,0)

eql[["Lorenzen"]]=
  lhEql(pars[["Lorenzen"]],range=c(min=0,max=120,minfbar=0,maxfbar=120,plusgroup=120))
eql[["0.2"]]     =
  lhEql(pars[["0.2"]],     range=c(min=0,max=120,minfbar=0,maxfbar=120,plusgroup=120))

ggplot(ldply(eql, function(x) as.data.frame(stock.wt(x),drop=T)))+
  geom_line(aes(age,data,col=.id))+
  theme_bw()+xlab("Age")+ylab("Kg")+
  theme(legend.position="none")

ggplot(ldply(eql, function(x) as.data.frame(m(x),drop=T)))+
  geom_line(aes(age,data,col=.id))

ggplot(ldply(eql, function(x) as.data.frame(stock.n(x)[,1],drop=T)))+
  geom_line(aes(age,data,group=.id))+
  facet_grid(.id~.,scale="free")

ggplot(ldply(eql, function(x) as.data.frame(stock.n(x)[,1]*stock.wt(x)[,1],drop=T)))+
  geom_line(aes(age,data,group=.id,col=.id))+
  facet_grid(.id~.,scale="free")+
  scale_x_continuous(limits=c(0,30))

aopt=ldply(eql,function(x) {
  ts=c(stock.n(x)[,1]*stock.wt(x))
  data.frame(age=rev(order(ts))[1]-1)})

lopt=ddply(aopt, .(.id), function(x) data.frame(length=c(FLife::vonB(FLQuant(x$age),par))))

lopt1=lopt(pars[["Lorenzen"]])
lopt2=lopt(pars[["0.2"]])

vonB(params=par,length=lopt1)
vonB(params=par,length=lopt2)

crv=mdply(data.frame(Length=seq(56,239)),function(Length){
  biomass=FLife::loptFn(Length,params=par)
  age    =c(vonB(length=FLQuant(Length),params=par))
  data.frame("Age"=age,"Biomass"=biomass)})

ggplot(crv)+
  geom_line(aes(Age,Biomass))
}


#res=params[1]%*%(a*(params["linf"]%*%(1.0-exp((-params["k"])%*%(age%-%params["t0"])))^b))%^%params[2])
