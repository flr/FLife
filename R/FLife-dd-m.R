if (FALSE){
  library(FLife)
  library(FLasher)
  
  lh=lhPar(FLPar(linf=100))
  eq=lhEql(lh)
  fbar(eq)[]=refpts(eq)["msy","harvest"]
  stk=as(eq,"FLStock")
  stk=fwd(stk,fbar=fbar(eq)[,-1]%=%refpts(eq)["msy","harvest"],sr=eq)
  
  fbar(eq)=FLQuant(c(sort(refpts(eq)[,"harvest"])),dimnames=list(year=1:5))
  nms=dimnames(sort(refpts(eq)[,"harvest"]))$refpt

  ggplot(exp(-lorenzenFn(stock.wt(eq))))+
    geom_line(aes(age,data,col=ac(year)))+
    xlab("Age")+ylab("Survival")

  dat=model.frame(FLQuants(survival=exp(-lorenzenFn(stock.wt(stk)[,1:5])),
                           n       =stock.n(eq)),drop=TRUE)
  
  surv=subset(dat,year==3)[,c("n","survival")]
  a=surv[,"survival"]
  v=subset(dat,year==1)[,c("n")]

  a=v[2]
  b=rep(.1,1)
  N=seq(1,100)
  
dd<-function(N,a,b){
  a*N/(1+N/b)}

ddm<-function(N,a,b){
  -log((a*N/(1+N/b))/N)}

cv<-function(x) var(x)^0.5/mean(x)

plot(dd(N,.1,25))

lh=lhPar(FLPar(linf=100))
eq=lhEql(lh)

srDev=rlnorm(1,FLQuant(0,dimnames=list(year=seq(1000))),0.3)

fbar(eq)=FLQuant(0,dimnames=list(year=seq(1000)))
fbar(eq)[]=refpts(eq)["msy","harvest"]*.1
stk=as(eq,"FLStock")
stk1=fwd(stk,fbar=fbar(eq)[,-1],sr=eq,residuals=srDev)

fbar(eq)[]=refpts(eq)["msy","harvest"]
stk=as(eq,"FLStock")
stk2=fwd(stk,fbar=fbar(eq)[,-1],sr=eq,residuals=srDev)

fbar(eq)[]=refpts(eq)["msy","harvest"]*2
stk=as(eq,"FLStock")
stk3=fwd(stk,fbar=fbar(eq)[,-1],sr=eq,residuals=srDev)

plot(FLStocks(list("0.2"=stk1,"msy"=stk2,"2"=stk3)))

stk1.=stk1
stk2.=stk2
stk3.=stk3

for (i in 2:1000){
  stk1@m[,i]=ddm(stock.n(stk1)[,i],surv[,"survival"],surv[,"n"])
  stk1=fwd(stk1,fbar=fbar(stk1)[,i],sr=eq,residuals=srDev)

  stk2@m[,i]=ddm(stock.n(stk2)[,i],surv[,"survival"],surv[,"n"])
  stk2=fwd(stk2,fbar=fbar(stk2)[,i],sr=eq,residuals=srDev)

  stk3@m[,i]=ddm(stock.n(stk3)[,i],surv[,"survival"],surv[,"n"])
  stk3=fwd(stk3,fbar=fbar(stk3)[,i],sr=eq,residuals=srDev)
  }

plot(FLStocks(list("0.2"=stk1[,ac(20:99)],
                   "msy"=stk2[,ac(20:99)],
                   "2"  =stk3[,ac(20:99)])))

x=rlnorm(100,0,0.3)
y10=dd(x*surv[2,"n"]*.5,surv[2,"survival"],surv[2,"n"])
y25=dd(x*surv[2,"n"],   surv[2,"survival"],surv[2,"n"])
y50=dd(x*surv[2,"n"]*2, surv[2,"survival"],surv[2,"n"])

ggplot(data.frame(x=rep(1:100,2),y=c(x,y),what=rep(c("input","output"),each=100)))+
  geom_line(aes(x,y,col=what))
}


#' mdd
#'
#' Lorenzen natural mortality relationship where M is a function of weight, modified to 
#' explicitly included M as a function of numbers in a cohort, i.e. density dependence
#'
#'  @details
#'    
#' The Lorenzen natural mortality relationship is a function of mass-at-age i.e. 
#'    M=a*wt^b 
#' 
#' The relationship can be explained by population density, since as fish grow they 
#' also die and so there is potentially less competition for resources between larger and 
#' older fish. Density dependence can be modelled by a logistic function, a sigmoid 
#' curve (or S shaped) curve, with equation 
#' 
#' f(x)=L/(1+exp(-k(x-x0)))
#'  
#' where
#' e  is the natural logarithm base (also known as Euler's number),
#' x0 is the x-value of the sigmoid's midpoint,
#' L  is the curve's maximum value, and
#' k  the steepness of the curve.
#' 
#' Combining the two functions gives
#'     
#'    M=aL/(1+exp(-k(n-ref)))*wt^b;
#'    
#' @param object  mass at which M is to be predicted
#' @param params an \code{FLPar} with two values; i.e. a equal to M at unit mass and b a power term; defaults are a=0.3 and b=-0.288
#' @param scale reference 
#' @param k rate of change in density dependence
#' @param m function with mortality model, by default gisalson
#' 
#' @aliases mdd mdd-method mdd,FLQuant,FLPar-method
#' @param ... other arguments, such as scale, e.g. stock numbers now relative to a reference level, e.g. at virgin biomass and k steepness of relationship
# 
#' @export
#' @docType methods
#' @rdname mdd
#' 
#' @seealso \code{\link{lorenzen}}
#'  
#' @examples
#' \dontrun{
#' library(FLBRP)
#' library(FLife)
#' 
#' data(teleost)
#' par=teleost[,"Hucho hucho"]
#' par=lhPar(par)
#' hutchen=lhEql(par)
#' 
#' scale=stock.n(hutchen)[,25]%*%stock.wt(hutchen)
#' scale=(stock.n(hutchen)%*%stock.wt(hutchen)%-%scale)%/%scale
#' 
#' m=mdd(wt2len(stock.wt(hutchen),par),params=par,scale,k=.9) 
#'  
#' ggplot(as.data.frame(m))+
#'    geom_line(aes(age,data,col=factor(year)))+
#'    theme(legend.position="none")+
#'    scale_x_continuous(limits=c(0,15))
#' 
#' m=mdd(stock.wt(hutchen),params=FLPar(m1=3,m2=-0.288),scale,k=1.2,m=lorenzen)   
#' }
setMethod('mdd', signature(object='FLQuant',params='FLPar'),
          function(object,params,scale,k=1,m=gislason) { 
            
            map=(2/(1+exp(k*scale)))
            
            res=map%*%m(object,params)
            res})

