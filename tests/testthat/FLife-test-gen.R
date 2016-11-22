library(FLife)

testObjects=list()

set.seed(123)
data(ple4)

#ages
testObjects[["ages"]]=ages(m(ple4))

#cc,numeric,numeric-method
ctc=as.data.frame(catch.n(ple4))
testObjects[["cc"]]=cc(age=ctc$age,n=ctc$data)


#dnormal
params=FLPar(a1=4,sl=2,sr=5000)
testObjects[["dnormal"]]=dnormal(FLQuant(1:10,dimnames=list(age=1:10)),params)


#gascuel
testObjects[["gascuel"]]=gascuel(10)

#gislason
params=lhPar(FLPar(linf=111))
len=FLQuant(c( 1.90, 4.23, 7.47,11.48,16.04,20.96,26.07,31.22,
               36.28,41.17,45.83,50.20,54.27,58.03,61.48,64.62),
             dimnames=list(age=1:16))
testObjects[["gislason"]]=gislason(len,params)


#gompertz
params=FLPar(linf=100,a=2,b=.4)
age=FLQuant(1:10,dimnames=list(age=1:10))
testObjects[["gompertz"]]=gompertz(age,params)

#knife
params=FLPar(a1=4)
age=FLQuant(1:10,dimnames=list(age=1:10))
testObjects[["knife"]]=knife(age,params)

#lambda
library(popbio)
eql=lhEql(lhPar(FLPar(linf=100)))
L=leslie(eql)
testObjects[["leslie"]]=lambda(L[drop=TRUE])



#len2wt
params=FLPar(a=1,b=3)
testObjects[["len2wt"]]=len2wt(FLQuant(10),params)


#lhEql
testObjects[["lhEql"]]=refpts(lhEql(lhPar(FLPar(linf=100))))["msy","yield"]


#lhPar
testObjects[["lhPar"]]=lhPar(FLPar(linf=200))

#lhRef
library(FLBRP)
params=FLPar(linf=100,t0=0,k=.4)
params=lhPar(params)
testObjects[["lhRef"]]=lhRef(params)

#logistic
params=FLPar(a50=4,ato95=1,asym=1.0)
age=FLQuant(1:10,dimnames=list(age=1:10))
testObjects[["logistic"]]=logistic(age,params)

#lopt
params=lhPar(FLPar(linf=100))
testObjects[["lopt"]]=lopt(params)

#loptAge
testObjects[["loptAge"]]=loptAge(params)

#lorenzen
mass=FLQuant(c( 1.90, 4.23, 7.47,11.48,16.04,20.96,26.07,31.22,
               36.28,41.17,45.83,50.20,54.27,58.03,61.48,64.62),
             dimnames=list(age=1:16))
testObjects[["lorenzen"]]=lorenzen(mass)


#matdd
testObjects[["matdd"]]=NULL
  

#mdd
eql=brp(FLBRP(ple4))
fbar(eql)=FLQuant(c(seq(0,4,length.out=101)*refpts(eql)["msy","harvest"]))
stk=as(eql,"FLStock")

fbar(eql)=fbar(eql)[,1] 
scale=(stock.n(stk)%-%stock.n(eql))%/%stock.n(eql)
par=FLPar(m1=.2,m2=-0.288)

testObjects[["m"]]=mdd(stock.wt(stk),par,scale)


x=rlnorm(100)
testObjects[["moment"]]=moment(x)

#noise
flq=FLQuant(1:100)
testObjects[["white"]]=noise(1000,flq,sd=.3,b=0)
testObjects[["red"]]  =noise(1000,flq,sd=.3,b=0.7)

#unbin
x=summary(cut(runif(100),seq(0,1,.1)))
testObjects[["unbin"]]=unbin(x)

#sigmoid
params=FLPar(a50=4,ato95=1)
age=FLQuant(1:10,dimnames=list(age=1:10))
testObjects[["sigmoid"]]=sigmoid(age,params)

#vonB
params=FLPar(linf=100,t0=0,k=.4)
age=FLQuant(1:10,dimnames=list(age=1:10))
len=vonB(age,params)
testObjects[[""]]=
#inverse growth curve
testObjects[["vonBInv"]]=vonB(params=params,length=len)

  
#wt2len
params=FLPar(a=0.1,b=3)
testObjects[["wt2len"]]=wt2len(FLQuant(10),params)

#sv
params=FLPar(a=37.8,b=8.93)
testObjects[["sv"]]=sv(params,"bevholt",.4)

#rod
object=rlnorm(1,FLQuant(0,dimnames=list(year=1:30)),.3)
rod=FLife:::rodFn
testObjects[["rod"]]=rod(object) 

