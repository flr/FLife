library(testthat)
library(FLife)

test_check("FLife")

### Name: ages,FLQuant-method
data(ple4)
ages(m(ple4))


### Name: cc,numeric,numeric-method
data(ple4)
ctc=as.data.frame(catch.n(ple4))
dat=cc(age=ctc$age,n=ctc$data)
head(dat)


#dnormal
params=FLPar(a1=4,sl=2,sr=5000)
dnormal(FLQuant(1:10,dimnames=list(age=1:10)),params)



#gascuel
gascuel(10)

#gislason
params=lhPar(FLPar(linf=111))
len=FLQuant(c( 1.90, 4.23, 7.47,11.48,16.04,20.96,26.07,31.22,
               36.28,41.17,45.83,50.20,54.27,58.03,61.48,64.62),
             dimnames=list(age=1:16))
gislason(length,params)



#gompertz
params=FLPar(linf=100,a=2,b=.4)
age=FLQuant(1:10,dimnames=list(age=1:10))
gompertz(age,params)

#knife
params=FLPar(a1=4)
age=FLQuant(1:10,dimnames=list(age=1:10))
knife(age,params)

#lambda
library(popbio)
eql=lhEql(lhPar(FLPar(linf=100)))
L=leslie(eql)
lambda(L[drop=TRUE])



#len2wt
params=FLPar(a=1,b=3)
len2wt(FLQuant(10),params)



#leslie
eql=lhEql(lhPar(FLPar(linf=100)))
leslie(eql)

#lhEql
eql=lhEql(lhPar(FLPar(linf=100)))


#lhPar
lhPar(FLPar(linf=200))

#lhRef
library(FLBRP)
params=FLPar(linf=100,t0=0,k=.4)
params=lhPar(params)
lhRef(params)

#logistic

params=FLPar(a50=4,ato95=1,asym=1.0)
age=FLQuant(1:10,dimnames=list(age=1:10))
logistic(age,params)




cleanEx()
nameEx("lopt")
#lopt

flush(stderr()); flush(stdout())

### Name: lopt,FLPar-method
### Title: Length at maximum biomass
### Aliases: lopt lopt,FLPar-method lopt-method

### ** Examples

## Not run: 
params=lhPar(FLPar(linf=100))
lopt(params)




cleanEx()
nameEx("loptAge")
#loptAge

flush(stderr()); flush(stdout())

### Name: loptAge,FLPar-method
### Title: Age at maximum biomass
### Aliases: loptAge loptAge,FLPar-method loptAge-method

### ** Examples

## Not run: 
params=lhPar(FLPar(linf=100))
loptAge(params)




cleanEx()
nameEx("lorenzen")
#lorenzen

flush(stderr()); flush(stdout())

### Name: lorenzen,FLQuant,FLPar-method
### Title: lorenzen
### Aliases: lorenzen lorenzen,FLQuant,FLPar-method
###   lorenzen,FLQuant,missing-method lorenzen,FLQuant,numeric-method
###   lorenzen,numeric,missing-method lorenzen-method

### ** Examples

## Not run: 
mass=FLQuant(c( 1.90, 4.23, 7.47,11.48,16.04,20.96,26.07,31.22,
               36.28,41.17,45.83,50.20,54.27,58.03,61.48,64.62),
             dimnames=list(age=1:16))
lorenzen(mass)




cleanEx()
nameEx("matdd")
#matdd

flush(stderr()); flush(stdout())

### Name: matdd,FLQuant,FLPar-method
### Title: matdd
### Aliases: matdd matdd,FLQuant,FLPar-method matdd-method

### ** Examples

## Not run: 
#bug

library(FLBRP)
library(FLife)

data(ple4)

sr=fmle(as.FLSR(ple4,model="bevholtSV"), fixed=list(s=.75,spr0=spr0(FLBRP(ple4))))
#control=list("silent"=TRUE))

eql=brp(FLBRP(ple4,sr=ab(sr)))
fbar(eql)=FLQuant(c(seq(0,4,length.out=101)*refpts(eql)["msy","harvest"]))
stk=as(eql,"FLStock")

fbar(eql)=fbar(eql)[,1]
scale=(stock.n(stk)%-%stock.n(eql))%/%stock.n(eql)
par=lhPar(FLPar(linf=100))

scale=exp(noise(1,m(stk),sd=.2,b=0.75))
 
ggplot(as.data.frame(mat,drop=TRUE))+
   geom_line(aes(age,data,group=year,col=factor(year)))+
   theme(legend.position="none")
   
 




cleanEx()
nameEx("mdd")
#mdd

flush(stderr()); flush(stdout())

### Name: mdd,FLQuant,FLPar-method
### Title: mdd
### Aliases: mdd mdd,FLQuant,FLPar-method mdd-method

### ** Examples

## Not run: 
library(FLBRP)
library(FLife)

data(ple4)

eql=brp(FLBRP(ple4))
fbar(eql)=FLQuant(c(seq(0,4,length.out=101)*refpts(eql)["msy","harvest"]))
stk=as(eql,"FLStock")

fbar(eql)=fbar(eql)[,1] 
scale=(stock.n(stk)%-%stock.n(eql))%/%stock.n(eql)
par=FLPar(m1=.2,m2=-0.288)

m=mdd(stock.wt(stk),par,scale)
ggplot(as.data.frame(m,drop=TRUE))+geom_line(aes(age,data,group=year,col=factor(year)))+theme(legend.position="none")


#moment
x=rlnorm(100)
moment(x)

#noise
flq=FLQuant(1:100)
white <- noise(1000,flq,sd=.3,b=0)

red <- noise(1000,flq,sd=.3,b=0.7)

#powh
x=summary(cut(runif(100),seq(0,1,.1)))
unbin(x)


#rod
 object=rlnorm(1,FLQuant(0,dimnames=list(year=1:30)),.3)
 pg=rod(object) 
 ggplot(object) +
     geom_polygon(aes(year,data,group=regime),
         fill="lavender",col="blue",
         lwd=.25,data=pg,alpha=.75)+
     geom_point(aes(year,data))+
     geom_line(aes(year,data))


#sigmoid
params=FLPar(a50=4,ato95=1)
age=FLQuant(1:10,dimnames=list(age=1:10))
sigmoid(age,params)

#sv

params=FLPar(a=37.8,b=8.93)
sv(params,"bevholt",.4)


#vonB
params=FLPar(linf=100,t0=0,k=.4)
age=FLQuant(1:10,dimnames=list(age=1:10))
len=vonB(age,params)

#inverse growth curve
vonB(params=params,length=len)

#wt2len
params=FLPar(a=0.1,b=3)
wt2len(FLQuant(10),params)


