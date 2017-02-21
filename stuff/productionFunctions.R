library(FLBRP)
library(FLife)
  
data(teleost)

alb=lhPar(teleost[,"Thunnus alalunga"])
bh9=lhEql(alb)

alb["s"]=0.9
bhDD=lhEql(alb)

fbarMult=c(rep(0.1,50),seq(0.1,4,length.out=100),seq(4,.7,length.out=10),rep(0.4,40))
fbar(bh9)=refpts(bh9)["msy","harvest"]*fbarMult
stk=fwd(bh9)

plot(stk)

fbar(bhDD)=refpts(bhDD)["msy","harvest"]*fbarMult
stkM=fwd(bhDD)
msyYr=sum(!(fbar(bhDD)>c(refpts(bhDD)["msy","harvest"])))
refYr=stock.n(bhDD)[,msyYr]%*%stock.wt(bhDD)

for (i in seq(dim(fbar(stkM))[2])[-1]){
  
  scale=(stock.n(stkM[,i])%*%stock.wt(stkM[,i])%-%refYr)%/%refYr
  
  m(stkM[,i])=mdd(wt2len(stock.wt(stkM[,i]),alb),params=alb,scale,k=0.1) 
  
  stkM=fwd(stkM,f=fbar(stkM[,i]),sr=bhDD)
  }

plot(FLStocks("M"=window(stkM,start=40),"SRR"=window(stk,start=40)))

alb[c("s","v")]=c(0.9,1000*max(catch(stk))/max(catch(stkM)))
bhDD=lhEql(alb)

fbar(bhDD)=refpts(bhDD)["msy","harvest"]*fbarMult
stkM=fwd(bhDD)
msyYr=sum(!(fbar(bhDD)>c(refpts(bhDD)["msy","harvest"])))
refYr=stock.n(bhDD)[,msyYr]%*%stock.wt(bhDD)

for (i in seq(dim(fbar(stkM))[2])[-1]){
  
  scale=(stock.n(stkM[,i])%*%stock.wt(stkM[,i])%-%refYr)%/%refYr
  
  m(stkM[,i])=mdd(wt2len(stock.wt(stkM[,i]),alb),params=alb,scale,k=0.1) 
  
  stkM=fwd(stkM,f=fbar(stkM[,i]),sr=bhDD)
}

#Maturity
alb=lhPar(teleost[,"Thunnus alalunga"])
alb["s"]=0.9
bhDD=lhEql(alb)

stkMat=fwd(bhDD)

fbar(bhDD)=refpts(bhDD)["msy","harvest"]*fbarMult
stkMat=fwd(bhDD)

msyYr=sum(!(fbar(bhDD)>c(refpts(bhDD)["msy","harvest"])))
refYr=stock.n(bhDD)[,msyYr]%*%stock.wt(bhDD)

for (i in seq(dim(fbar(stkMat))[2])[-1]){
  
  scale=(stock.n(stkMat[,i])%*%stock.wt(stkMat[,i])%-%refYr)%/%refYr
  
  mat(stkMat[,i])[]=matdd(ages(stock.wt(stkMat[,i])),params=alb,scale,k=2.0) 
  
  stkMat=fwd(stkMat,f=fbar(stkMat[,i]),sr=bhDD)
}

plot(FLStocks("Mat"=window(stkMat,start=40),"SRR"=window(stk,start=40)))

alb=lhPar(teleost[,"Thunnus alalunga"])
alb[c("s","v")]=c(0.9,1000*max(catch(stk))/max(catch(stkMat)))
bhDD=lhEql(alb)

fbar(bhDD)=refpts(bhDD)["msy","harvest"]*fbarMult
stkMat=fwd(bhDD)

msyYr=sum(!(fbar(bhDD)>c(refpts(bhDD)["msy","harvest"])))
refYr=stock.n(bhDD)[,msyYr]%*%stock.wt(bhDD)

for (i in seq(dim(fbar(stkMat))[2])[-1]){
  
  scale=(stock.n(stkMat[,i])%*%stock.wt(stkMat[,i])%-%refYr)%/%refYr
  
  mat(stkMat[,i])[]=matdd(ages(stock.wt(stkMat[,i])),params=alb,scale,k=2.0) 
  
  stkMat=fwd(stkMat,f=fbar(stkMat[,i]),sr=bhDD)
}

plot(FLStocks("M"  =window(stkM,  start=40),
              "Mat"=window(stkMat,start=40),
              "SRR"=window(stk,   start=40)))

dat=as.data.frame(FLQuants("M"  =catch(stkM)/fbar(stkM),
                           "Mat"=catch(stkMat)/fbar(stkMat),
                           "SRR"=(catch(stk)/fbar(stk))))
dat=cbind(dat,fmsy=rep(c(rep(0.1,50),seq(0.1,4,length.out=100)),3))

ggplot(subset(dat,fmsy>0.1&fmsy<2))+
  geom_line(aes(fmsy,log(data),col=qname))+
  xlab(expression(F[MSY]))+ylab("CPUE")

library(ggplotFL)
library(FLAssess)
stk.=stkM
m(      stk.)=m(stk)
stk.=stk.+VPA(stk.,fratio=1)
plot(FLStocks(M=stkM,SRR=stk.))    

sep.vpa.control <- FLSepVPA.control(sep.nyr=3, sep.age=2, sep.sel=.7)

stk..<-SepVPA(stk.,sep.vpa.control,fratio=1)
harvest(stk.)=harvest(stk..)
stock.n(stk.)=stock.n(stk..)

plot(FLStocks(M=stkM,SRR=stk.))    



