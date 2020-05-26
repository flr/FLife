### Powers/Forrest
library(FLCore)
library(FLasher)
library(FLBRP)
library(FLXSA)
library(ggplotFL)
library(FLife)
library(plyr)

## Calculate numbers and catch-at-age given density depemdence
ddNFunc<-function(object,A=0.1){
  
  res=window(stock.n(object),end=dim(stock.n(object))[2]+1)
  res=setPlusGroup(res,dims(stock.n(object))$max+1)
  res[-1,-1]=NA
  
  a= stock.n(object)[,1]*exp(-object@m[,1]-harvest(object)[,1])
  b=1+stock.n(object)[,1]*(A)/(object@m[,1]+harvest(object)[,1])*(1-exp(-object@m[,1]-harvest(object)[,1]))
  
  res[-1,2]=a/b
  
  for (i in seq(dim(stock.n(object))[2])){
    
    ## N+1
    a= res[-dim(res)[1],i]*exp(-object@m[,i]-harvest(object)[,i])
    b=1+res[-dim(res)[1],i]*(A)/(object@m[,i]+harvest(object)[,i])*(1-exp(-object@m[,i]-harvest(object)[,i]))
  
    res[-1,i+1]=a/b}
  
  object=object[,-1]
  res   =res[-1,-c(1,dim(res)[2])]
  
  stock.n(object)[-1,]=res[-dim(res)[1]]
  stock.n(object)[dim(res)[1],]=res[dim(res)[1],]+res[dim(res)[1],]
  
  ## Catch
  catch.n(object)=
    harvest(object)/A*
    log(1+stock.n(object)*A/(object@m+harvest(object))*(1-exp(-object@m-harvest(object))))
  
  catch(object)=computeCatch(object)
  
  return(object)}

# M given demsity dependence
mFunc<-function(object,A){
  if (is.numeric(A)) A=FLPar(A)
  
  object@m%+%(A%*%stock.n(object))}  

# M from N-at-age
mFunc2<-function(object){
  res=object@m
  res[-dim(stock.n(object))[1],-dim(stock.n(object))[2]]=
    log(stock.n(object)[-dim(stock.n(object))[1],-1]/
          stock.n(object)[-1,-dim(stock.n(object))[2]])-
    harvest(object)[-dim(stock.n(object))[1],-dim(stock.n(object))[2]]
  res[dim(stock.n(object))[1],]=res[dim(stock.n(object))[1]-1,]
  res[,dim(stock.n(object))[2]]=res[,dim(stock.n(object))[1]-1]
  
  res}  

# Run XSA
xsaMP<-function(om,pg=10,ctrl=FLXSA.control(maxit=250)){
  stk  =setPlusGroup(om,pg)
  idx  =FLIndex(index=stock.n(stk))
  range(idx)[c("plusgroup","startf","endf")]=c(pg,0.1,.2)
  stk+FLXSA(stk,idx,control=ctrl,diag.flag=FALSE)}

# Create FLstock with constant M

mVector<-function(linf,mval,A){
  lh=lhPar(FLPar(linf=linf),s=1,m=list(model="constant",params=c(m1=mval)))
  eq=lhEql(lh)
  fbar(eq)=FLQuant(c(refpts(eq)["msy","harvest"]))
  
  A =FLQuant(A,dimnames=dimnames(eq@m))
  
  a= stock.n(eq)[,1]*exp(-eq@m[,1]-harvest(eq)[,1])
  b=1+stock.n(eq)[,1]*(A)/(eq@m[,1]+harvest(eq)[,1])*(1-exp(-eq@m[,1]-harvest(eq)[,1]))
  
  skn=stock.n(eq)
  sdd=a/b
  skn[dim(skn)[1]-1]=sdd[dim(skn)[1]]+sdd[dim(skn)[1]-1]
  skn[-1]=sdd[dim(skn)[1]-1]
  
  eq@m%+%(A*skn)}

## stock with declining M-at-age
stkM  =lhStk(FLPar(linf=50),
           s    =1,
           m    =mVector(50,0.2,1/50),
           fmult=function(x) FLQuant(rep(1,1001))%*%refpts(x)["f0.1","harvest"])

# Recruitment noise
srDev    =rlnorm(1,rec(stkM[[1]])%=%0,0.5)
stkM[[1]]=fwd(stkM[[1]],fbar=fbar(stkM[[2]])[,-1],sr=stkM[[2]],residuals=srDev)

stk=stkM[[1]]
stk@m[]=0.2
stkdd=ddNFunc(stk,1/50)

plot(FLStocks("M"=stkM[[1]],"DD"=stkdd))
plot(FLQuants("M" =ssb(stkM[[1]])[,ac(201:300)],
              "DD"=ssb(stkdd)[   ,ac(201:300)]))

stkMP=stkM[[1]][,ac(2:1001)]
catch.n(stkMP)=catch.n(stkdd)[,ac(2:1001)]
catch(stkMP)=computeCatch(stkMP)
stkMP=xsaMP(stkMP,pg=10)
 
plot(FLStocks("M"=stkM[[1]],"DD"=stkdd,"MP"=stkMP))
plot(FLQuants("M" =ssb(stkM[[1]])[,ac(201:300)],
              "DD"=ssb(stkdd)[    ,ac(201:300)],
              "MP"=ssb(stkMP)[    ,ac(201:300)]))

plot(FLQuants("DD"=rec(stkdd)[    ,ac(201:300)],
              "MP"=rec(stkMP)[    ,ac(201:300)]))

## Stock with DD M
#stk=lhStk(FLPar(linf=50),s=1,m=list(model="constant",params=c(m1=0.2)),
#            fmult=function(x) FLQuant(c(refpts(x)["msy","harvest"]),dimnames=list(year=seq(1001))))[[1]]


############# Verhurlst
  lh=lhPar(FLPar(linf=20))
  eq=lhEql(lh)
  fbar(eq)[]=refpts(eq)["f0.1","harvest"]
  stk=as(eq,"FLStock")
  stk=fwd(stk,fbar=fbar(eq)[,-1]%=%refpts(eq)["f0.1","harvest"],sr=eq)
  
  fbar(eq)=FLQuant(c(sort(refpts(eq)[,"harvest"])),dimnames=list(year=1:5))
  nms=dimnames(sort(refpts(eq)[,"harvest"]))$refpt
  
  ggplot(exp(-lorenzenFn(stock.wt(eq))))+
    geom_line(aes(age,data,col=ac(year)))+
    xlab("Age")+ylab("Survival")
  
  dat=model.frame(FLQuants(survival=exp(-lorenzenFn(stock.wt(stk)[,1:5])),
                           n       =stock.n(eq)),drop=TRUE)
  
  surv=subset(dat,year==3)[,c("n","survival")]
  
  dd<-function(N,a,b){
    a*N/(1+N/b)}
  
  ddm<-function(N,a,b){
    -log((a*N/(1+N/b))/N)}
  
  cv<-function(x) var(x)^0.5/mean(x)
  
  plot(ddm(N,.1,50))
  
  lh=lhPar(FLPar(linf=20))
  eq=lhEql(lh)
  
  srDev=rlnorm(1,FLQuant(0,dimnames=list(year=seq(1000))),0.3)
  
  fbar(eq)=FLQuant(0,dimnames=list(year=seq(1000)))
  fbar(eq)[]=refpts(eq)["f0.1","harvest"]*.1
  stk=as(eq,"FLStock")
  stk1=fwd(stk,fbar=fbar(eq)[,-1],sr=eq,residuals=srDev)
  
  fbar(eq)[]=refpts(eq)["f0.1","harvest"]
  stk=as(eq,"FLStock")
  stk2=fwd(stk,fbar=fbar(eq)[,-1],sr=eq,residuals=srDev)
  
  fbar(eq)[]=refpts(eq)["f0.1","harvest"]*2
  stk=as(eq,"FLStock")
  stk3=fwd(stk,fbar=fbar(eq)[,-1],sr=eq,residuals=srDev)
  
  plot(FLStocks(list("0.2"=stk1,"f0.1"=stk2,"2"=stk3)))
  
  stk1.=stk1
  stk2.=stk2
  stk3.=stk3
  
  for (i in 101:999){
    stk1@m[,i]=ddm(stock.n(stk1)[,i],surv[,"survival"],surv[,"n"]*2.5)
    stk1=fwd(stk1,fbar=fbar(stk1)[,i],sr=eq,residuals=srDev)
    
    stk2@m[,i]=ddm(stock.n(stk2)[,i],surv[,"survival"],surv[,"n"]*2.5)
    stk2=fwd(stk2,fbar=fbar(stk2)[,i],sr=eq,residuals=srDev)
    
    stk3@m[,i]=ddm(stock.n(stk3)[,i],surv[,"survival"],surv[,"n"]*2.5)
    stk3=fwd(stk3,fbar=fbar(stk3)[,i],sr=eq,residuals=srDev)
  }
  
  plot(FLStocks(list("0.2"=stk1[,ac(20:999)],
                     "f0.1"=stk2[,ac(20:999)],
                     "2"  =stk3[,ac(20:60)])))
  
  fb=fbar(eq)
  for (i in 1:100){
    fbar(eq)=FLQuant(c(fb[,i]))
    
    eq@m[]=ddm(stock.n(eq)[,i],surv[,"survival"],surv[,"n"])
    m
  }
  
  x=rlnorm(100,0,0.3)
  y10=dd(x*surv[2,"n"]*.5,surv[2,"survival"],surv[2,"n"])
  y25=dd(x*surv[2,"n"],   surv[2,"survival"],surv[2,"n"])
  y50=dd(x*surv[2,"n"]*2, surv[2,"survival"],surv[2,"n"])
  
  ggplot(data.frame(x=rep(1:100,2),y=c(x,y),what=rep(c("input","output"),each=100)))+
    geom_line(aes(x,y,col=what))
}


