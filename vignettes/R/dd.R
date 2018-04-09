library(FLife)
library(FLRP)

par=FLPar(teleost[c("linf","k","t0","l50","a","b"),"Thunnus alalunga"])
par=lhPar(par)
alb=lhEql(par)

fbar(alb)=FLQuant(c(seq(0,4,length.out=101)))*refpts(alb)["msy","harvest"]
names(dimnames(fbar(alb)))[1]="age"
stk=as(alb,"FLStock")
stk=fwd(stk,f=fbar(alb)[,-1],sr=alb)

ref=stock.n(alb)[,26]
ms=mdply(data.frame(i=seq(dims(fbar(alb))$year)[-1]),
          function(i){
              scale=((stock.n(stk)[,i]%-%ref)%/%ref)^0.1
              #length=exp(log(stock.wt(alb)%/%par["a"])%/%par["b"])
              as.data.frame(mdd(stock.wt(alb)/200,par,scale,k=.7))
              })

ggplot(as.data.frame(m(stk)))+
  geom_line(aes(age,data,col=factor(year)))+
  theme(legend.position="none")+
  scale_x_continuous(limits=c(0,15))


p=ggplot(mDD)+geom_line(aes(ssb,catch))

y  =eql[["Cushing"]]
ref=stock.n(y)[,1]
n  =stock.n(y)[,1]
mat=mat(y)#[,max(i-1,dims(m(y))$year)]    

matDD=mdply(data.frame(i=seq(dims(fbar(eql[["Cushing"]]))$year)[-1]),
            function(i){
              
              repeat{
                scale=(stock.n(y)[,i]%-%ref)%/%ref
                mat(y)[] =matdd(ages(scale),alb,scale,k=.9)
                
                if(sum((mat-mat(y))^2)<1e-6){break}
                mat=mat(y)}
              
              subset(model.frame(FLQuants(y,"catch","ssb"),drop=TRUE),
                     year==i)})
```


```{r fig6,fig.height=5,fig.width=6}
data=rbind.fill(
  cbind("Form"="M",mDD),cbind("Form"="Fecundity",matDD),
  cbind("Form"="Stock-recruit",model.frame(FLQuants(eql[["Cushing"]],"catch","ssb"),drop=TRUE)))

ggplot(data)+
  geom_line(aes(ssb,catch,colour=Form))+theme_bw()+theme(legend.position="bottom")+
  scale_colour_manual("",values=c("red","green","blue"))+
  xlab("Biomass")+ylab("Yield")
```

```{r fig7,fig.height=6}

k   =1
sr  ="cushing"
s   =0.7
alb  =FLPar(unlist(teleost[teleost$species=="Thunnus alalunga",
                           c("linf","k","t0","l50","a","b")]))
alb  =lhPar(rbind(alb,FLPar(m1=0.3,m2=-0.288,s=s)))
eq      =lhEql(alb,m=lorenzen,sr=sr)
fbar(eq)=FLQuant(c(seq(0,5,length.out=51)*
                     refpts(eq)["msy","harvest"]))
f       =fbar(eq)[,-1]

stkr=as(eq,"FLStock")
stkr=fwd(stkr,f=f,sr=eq)
ref =stock.n(eq)[,1]

stkm=stkr
m(  stkm)=mdd(stock.wt(stkm),alb,(stock.n(stkm)%-%ref)%/%ref,k) 
stkm=fwd(stkm,f=f,sr=eq)

stkf=stkr
mat(stkf)=matdd(ages(stock.wt(stkf)),alb,
                (stock.n(stkf)%-%ref)%/%ref,k,TRUE) 
stkf=fwd(stkf,f=f,sr=eq)

for (i in seq(dims(fbar(eq))$year)[2:50]){
  scale =(stock.n(stkm)[,i]%-%ref)%/%ref
  m(stkm)[,i]=mdd(stock.wt(stkm)[,i],alb,scale,k)
  
  stkm=fwd(stkm,f=f[,i],sr=eq)}

for (i in seq(dims(fbar(eq))$year)[2:50]){
  scale =(stock.n(stkf)[,i]%-%ref)%/%ref
  mat(stkf)[,i]=matdd(ages(scale),alb,scale,k,TRUE)
  
  stkf=fwd(stkf,f=f[,i],sr=eq)}

plot(FLStocks("SRR"=stkr[,-(1:1)],
              "M"  =stkm[,-(1:1)],
              "Fecundity"=stkf[,-(1:1)]))+
  theme_bw()+
  theme(legend.position="bottom")

p=ggplot(as.data.frame(FLQuants("SRR"=ssb(stkr),
                                "M"=ssb(stkm),
                                "Fecundity"=ssb(stkf))))+
  geom_line(aes(year,data,col=qname))+
  theme_bw()+
  theme(legend.position="bottom")
```


```{r fig8,fig.height=6}
k   =1
sr  ="cushing"
s   =0.7
alb  =FLPar(unlist(teleost[teleost$species=="Thunnus alalunga",
                           c("linf","k","t0","l50","a","b")]))
alb  =lhPar(rbind(alb,FLPar(m1=0.3,m2=-0.288,s=s)))
eq      =lhEql(alb,m=lorenzen,sr=sr)
fbar(eq)=FLQuant(c(seq(5,.5,length.out=51)*
                     refpts(eq)["msy","harvest"]))
f       =fbar(eq)[,-1]

stkr=as(eq,"FLStock")
stkr=fwd(stkr,f=f,sr=eq)
ref =stock.n(eq)[,1]

stkm=stkr
m(  stkm)=mdd(stock.wt(stkm),alb,(stock.n(stkm)%-%ref)%/%ref,k) 
stkm=fwd(stkm,f=f,sr=eq)

stkf=stkr
mat(stkf)=matdd(ages(stock.wt(stkf)),alb,
                (stock.n(stkf)%-%ref)%/%ref,k,TRUE) 
stkf=fwd(stkf,f=f,sr=eq)

for (i in seq(dims(fbar(eq))$year)[2:50]){
  scale =(stock.n(stkm)[,i]%-%ref)%/%ref
  m(stkm)[,i]=mdd(stock.wt(stkm)[,i],alb,scale,k)
  
  stkm=fwd(stkm,f=f[,i],sr=eq)}

for (i in seq(dims(fbar(eq))$year)[2:50]){
  scale =(stock.n(stkf)[,i]%-%ref)%/%ref
  mat(stkf)[,i]=matdd(ages(scale),alb,scale,k,TRUE)
  
  stkf=fwd(stkf,f=f[,i],sr=eq)}

plot(FLStocks("SRR"=stkr[,-(1:1)],
              "M"  =stkm[,-(1:1)],
              "Fecundity"=stkf[,-(1:1)]))+
  theme_bw()+
  theme(legend.position="bottom")

p=ggplot(as.data.frame(FLQuants("SRR"=ssb(stkr),
                                "M"=ssb(stkm),
                                "Fecundity"=ssb(stkf))))+
  geom_line(aes(year,data,col=qname),size=2)+
  theme(legend.position="bottom")
```
