power<-function(freq=1,cv=seq(0.2,0.4,0.1),nyr=30,pr=500,trend=2){
  
  #require(fishmethods)
  #require(akima)
  
  pwr=mdply(expand.grid(pse=cv,yrs=seq(nyr)), function(pse,yrs) {
    t.=powertrend(trend=trend,A1=1,PSE=pse,maxyrs=yrs,step=1,graph=FALSE,pR=pr)
    
    subset(t.[,c("R","power")], !(is.na(power)))})
  
  ### Power Analysis contours
  pws=mdply(data.frame(freq=freq),function(freq) transform(pwr, r=log(1+R/100)/yrs))
  pws=transform(subset(pws, power<.7 & power>.5), increase=ifelse(r>0,TRUE,FALSE))
  
  p.=ddply(subset(pws,!(is.na(r)|is.infinite(r))),.(freq,pse,increase), function(dat,power=0.6){  
    # linear interpolation at points (1,2), (5,6) and (10,12)
    akima.lip<-interpp(dat$yrs, dat$power, dat$r, seq(1,nyr,1),rep(power,nyr),duplicate="mean")
    data.frame(Year=akima.lip$x,r=akima.lip$z)})
  
  transform(p.,grp=paste(freq,pse,ac(increase)),cv=factor(pse),interval=factor(freq))}

if (FALSE){
  library(ggplot2)
  library(plyr)
  library(FLCore)
  library(ggplotFL)
  library(FLRP)
  library(FLash)
  
  
dat=power()

ggplot()+
  geom_line(aes(Year,r,group=grp,col=cv),data=dat)+
  geom_hline(aes(yintercept=seq(-.5,.25,.25)),linetype=2) + 
  geom_vline(aes(xintercept=seq(1,30,1)),linetype=21)+
  ylab("Population Growth Rate (r)")        +
  xlab("Year of Survey")

load("/home/laurie/Desktop/ownCloud-2016/MedSWOAss/MedSWOAss/Analysis/xsa/data/swom.RData")
load("/home/laurie/Desktop/ownCloud-2016/MedSWOAss/MedSWOAss/Analysis/xsa/data/eql.RData")

swo   =fwdWindow(window(swom,end=2013),eql[[1]],end=2045)
target=FLQuants(mlply(data.frame(F=seq(0,0.5,0.05)*refpts(eql[[1]])["msy","harvest"]), 
                      function(F) FLQuant(F,dimnames=list(year=2013:2045))))
names(target)=rep("f",11)

swo=fwd(swo,target,sr=eql[[1]])

plot(swo)

r=ldply(swo,function(x) log(c(ssb(x)[,"2017"]%/%ssb(x)[,"2012"]))/5)
r=transform(r,F=(seq(0,0.5,0.05)*refpts(eql[[1]])["msy","harvest"])[X1])[,-1]
names(r)=c("r","F")

dat=cbind(dat,F=predict(lm(F~r,data=r),data.frame(r=dat$r)))

p=ggplot()+
  geom_line(aes(Year,r,group=grp,linetype=cv),data=dat)+
  geom_hline(aes(yintercept=r),col="grey80",data=data.frame(r=seq(-0.4,0.4,0.1)))+ 
  geom_vline(aes(xintercept=seq(0,30,5)),linetype=1,col="grey80",size=.5)+
  ylab("Population Growth Rate (r)")        +
  xlab("Years of Survey")+
  theme_bw(14)+
    scale_x_continuous(breaks=seq(0,30,5))+
    scale_y_continuous(limits=c(-0.32,0.32),breaks=seq(-0.3,.3,.1),
                       labels=c("-0.3","-0.2","-0.1","0.0","0.1","0.2","0,3"))

ggsave("/home/laurie/Desktop/powerFig.jpeg",device="jpeg",dpi=1000)  
predict(lm(F~r,data=dat),new=data.frame(r=seq(-0.2,0.2,0.1)))
}