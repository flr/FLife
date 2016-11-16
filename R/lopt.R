globalVariables(c("aaply"))


#' lopt
#'
#' @description Finds length at maximum biomass, assuminmg natural mortality is a function of
#'  mass-at-age, i.e. Lorenzen. 
#' 
#' 
#' @param params FLPar
#' @param m A function, i.e. gislason
#' @param ... any other arguments
#' 
#' @aliases lopt lopt-method lopt,FLPar-method
#' 
#' @return FLPar with length at maximum biomass 
#' 
#' @details There are several ways to calculate \deqn{L_{opt}}, i.e.
#' i) \deqn{{2/3}^{rds}  L_{\infty}}
#' ii) \deqn{L_{\infty}\frac{3}{3+k/m}}
#' iii) by maximising the biomass of
#' iv) from an FLBRP object by fishing at F=0 and finding age where biomass is a maximum
#' 
#' @export
#' @docType methods
#' @rdname lopt
#' 
#' @seealso \code{\link{vonB}}  
#' 
#' @examples
#' \dontrun{
#' data(pars)
#' lopt(pars[[1]])
#' }
setMethod("lopt", signature(params="FLPar"),
       function(params,
                m=function(length,params) exp(0.55)*(length^-1.61)%*%(params["linf"]^1.44)%*%params["k"],
                ...){   
            dmns=dimnames(params)
            dmns$params="lopt"
            dm  =dim(params)
            
            res=aaply(params,seq(length(dm))[-1],function(x){
                   x.=FLPar(x)

                   rtn=try(optimise(loptFn,c(.01,c(x["linf"])*.99),params=x.,maximum=TRUE,m=m)$maximum)
   
                   if ("character" %in% mode(rtn)) rtn=NA
                   rtn})
            
            FLPar(array(res,dim=c(1, dm[-1]),dimnames=dmns))})

loptFn=function(x,params,
                m=function(length,params) exp(0.55)*(length^-1.61)%*%(params["linf"]^1.44)%*%params["k"],
                age=0:200){
  
 
  length=vonB(FLQuant(age,dimnames=list(age=age)),params)
  m.    =FLQuant(m(length,params), dimnames=list(age=age))
  mCum  =FLQuant(aaply(m.,2,cumsum),dimnames=dimnames(m.))

  a =qmax(vonB(params=params,length=FLQuant(x)),min(age))
  a =qmin(a,max(age))
  
  aMin=floor(a)
  if (is.na(aMin)) return (as.numeric(NA))
  
  aMax=ceiling(a)
  if (is.na(aMax)) return (as.numeric(NA))

  m1=m(vonB(aMin,params),params)
  m2=m(vonB(aMax,params),params)

  slope    =(m1-m2)/(aMin-aMax)
  intercept=m1-slope*aMin
  m.=(intercept+slope*a)*(a-aMin)
    
  n =exp(-mCum[ac(aMin)]-m.)
  
  c(n%*%len2wt(x,params))
  }

#params=par[,10]
#params["t0"]=-.1
#loptFn(x,params)

loptFn2<-function(params) params["linf"]*3/(3+exp(params["m2"])/params["k"])

#' loptAge
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
#' @return FLPar with length at maximum biomass 
#' 
#' @details There are several ways to calculate \deqn{L_{opt}}, i.e.
#' i) \deqn{{2/3}^{rds}  L_{\infty}}
#' ii) \deqn{L_{\infty}\frac{3}{3+k/m}}
#' iii) by maximising the biomass of
#' iv) from an FLBRP object by fishing at F=0 and finding age where biomass is a maximum
#' 
#' @export
#' @docType methods
#' @rdname loptAge
#' 
#' @seealso \code{\link{vonB}}  
#' 
#' @examples
#' \dontrun{
#' data(pars)
#' lopt(pars[[1]])
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
        m.    =FLQuant(m(length,params),    dimnames=dmns)
        mCum  =FLQuant(aaply(m.,6,sum))
        n     =exp(-mCum)
        c(n*FLife::len2wt(length[ac(ceiling(x))],params))}
           
      dmns=dimnames(params)
      dmns$params="lopt"
      dm  =dim(params)
        
      res=aaply(params,seq(length(dm))[-1],function(x){
            x.=FLPar(x)
            rtn=try(optimise(loptFn,c(0,40),params=x.,maximum=TRUE,m=m)$maximum)
            if ("character" %in% mode(rtn)) rtn=NA
            rtn})
     
      vonB(as(FLPar(array(res,dim=c(1, dm[-1]),dimnames=dmns)),"FLQuant"),params)
      })

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
