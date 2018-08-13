utils::globalVariables(c("rmvnorm"))
utils::globalVariables(c("bin"))
utils::globalVariables(c("V1"))
utils::globalVariables(c("bin"))
utils::globalVariables(c("p"))
utils::globalVariables(c("mp"))
utils::globalVariables(c("mp"))
utils::globalVariables(c("powertrend"))
utils::globalVariables(c("R"))
utils::globalVariables(c("yrs"))
utils::globalVariables(c("pse"))
utils::globalVariables(c("increase"))
utils::globalVariables(c("interpp"))
utils::globalVariables(c("regime"))
utils::globalVariables(c("b.msy"))
utils::globalVariables(c("b.current"))
utils::globalVariables(c("spectraFLQ"))

fn=function(year,object,sigma,nsample=100) {
  n  =as.data.frame(stock.n(object)[,year]^(1/3),drop=TRUE)
  len=c((stock.wt(object)[,year])^(1/3))
  len=exp(rmvnorm(n=nsample, mean=log(len), sigma=sigma))
  len=melt(len)
  names(len)=c("iter","age","length")
  res=merge(len,n,by="age")
  
  lmat=c(min(stock.wt(object)[,1][mat(object)[,1]>=0.5])^(1/3))
  
  idx =transform(res, bin=cut(length,c(0,lmat,lmat*1.5,1000)))
  idx =ddply(idx,.(bin),with,sum(data))
  idx =transform(idx,p=V1/sum(V1))
  idx}

utils::globalVariables(c("adply"))

fn2<-function(object,cv=0.3,nsample=100){
  
  sigma=matrix(0, ncol=dim(m(object))[1],nrow=dim(m(object))[1])
  diag(sigma)=cv

  mdply(data.frame(year=dimnames(m(object))$year),fn,object=object,sigma=sigma,nsample=nsample)}
  
ldex<-function(object,cv=0.3,nsample=100){
  
  idx=mdply(data.frame(iter=seq(dim(stk)[6])), 
            function(iter) fn2(iter(object,iter),cv=cv,nsample=nsample))
  idx=FLQuants(dlply(idx,.(bin),with, 
            as.FLQuant(data.frame(year=as.numeric(ac(year)),iter=iter,data=p))))
  
  idx[[2]]=idx[[2]]+idx[[3]]
  names(idx)=c("immature","mature","mega")}

if (FALSE){
  library(ggplot2)
  library(ggplotFL)
  library(FLBRP)
  library(FLife)
  library(mvtnorm)
  
par=lhPar(FLPar(linf=100))
stk=fwd(lhEql(par))[,ac(10:90)]
stk=setPlusGroup(stk,8)[-1,]
stk=propagate(stk,10)

plot(stk)
idx=ldex(stk)

ggplot(as.data.frame(idx[[1]]))+
  geom_boxplot(aes(as.factor(year),data))

ggplot(as.data.frame(idx[[2]]+idx[[2]]))+
  geom_boxplot(aes(as.factor(year),data))

ggplot(as.data.frame(idx[[3]]))+
  geom_boxplot(aes(as.factor(year),data))
}

if (FALSE){
  library(FLBRP)
  library(FLife)
  
  par=lhPar(FLPar(linf=100))
  stk=fwd(lhEql(par))
  stk=setPlusGroup(stk,8)[-1,]
  
  lo=c(lopt(par))
  lm=c(vonB(par,age=FLQuant(par["a50"])))
  
  sigma=matrix(0, ncol=8,nrow=8)
  diag(sigma)=0.3
  
  sim=mdply(data.frame(year=1:100), function(year) {
    n  =as.data.frame(stock.n( stk)[,year]^(1/3),drop=TRUE)
    len=c((stock.wt(stk)[,year]/par["a"])^(1/3))
    len=exp(rmvnorm(n=102, mean=log(len), sigma=sigma))
    len=melt(len)
    names(len)=c("iter","age","length")
    merge(len,n,by="age")})
  
  ggplot(mutate(sim,decade=10*((year-1)%/%10),yr=year-decade-1))+
    geom_histogram(aes(length,weight=data))+
    geom_vline(aes(xintercept=lm),col="blue")+
    geom_vline(aes(xintercept=lo),col="red")+
    geom_vline(aes(xintercept=lm*2),col="green")+
    facet_wrap(~decade)+
    scale_x_continuous(limits=c(0,125))
  
  sim=transform(sim, bin=cut(sim$length,c(0,lm,lo,lm*1.5,1000)))
  idx=ddply(sim,.(year,bin),with,sum(data))
  idx=ddply(idx,.(year),transform,p=V1/sum(V1))
  ggplot(aes(year,p,col=bin),data=idx)+
    geom_point()+
    geom_smooth()+
    facet_grid(bin~.,scale="free")
}

