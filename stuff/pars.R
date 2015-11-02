library(ggplot2)
library(FLCore)
library(FLBRP)
library(ggplotFL)

source('~/Desktop/flr/git/FLife/R/lh-growth.R')
source('~/Desktop/flr/git/FLife/R/gascuel.R')
source('~/Desktop/flr/git/FLife/R/lh-m-lorenzen.R')
source('~/Desktop/flr/git/FLife/R/lh-growth-vonB.R')
source("/home/laurie/Desktop/flr/git/FLife/R/lh-func.R")
source("/home/laurie/Desktop/flr/git/FLife/R/lh-gislasim.R")
source("/home/laurie/Desktop/flr/git/FLife/R/lh-ogive-dnormal.R")

##Atlantic Bigeye
betPar=FLPar(c(linf=217.3,   k  =0.18,  t0=-0.709,
               a   =2.396e-5,b  =2.9774,
               a50 =3,       l50=105,
               m1  =0.2795,  m2 =-0.288))

##Atlantic Yellowfin
#M assumed to be 0.8 for ages 0 and 1, and 0.6 for ages 2+
#Birthday February 14
#PG Age 5+ 
#Mat knife-edge at 3
yftPar=FLPar(c(linf=137.0,k=0.808,t0=-0.01,
               a=2.1527e-5,b=2.976,
               a50=3,l50=NA,
               gas.a=37.8,gas.b=8.93,gas.c=137.0,
               gas.d=8.93,gas.e=0.808,gas.f=7.49,
               m1   =0.3,m2=-0.288))
yftPar["l50"]=gascuel(yftPar,yftPar["a50"])

m =FLQuant(c(rep(0.8,2),rep(0.6,11)),dimnames=list(age=0:12))
ln=gascuel(yftPar,ages(m))
wt=len2wt(yftPar,ln)
m.=lorenzen(yftPar[c("m1","m2")],wt)

fn=function(m1,m,par) {
  par["m1"]=m1
  m.=lorenzen(par,wt)  
  sum((m-m.)^2,na.rm=T)}

yftPar["m1"]=optimise(fn, c(0.001,10),m=m,par=yftPar[c("m1","m2")])$minimum

m.=lorenzen(yftPar[c("m1","m2")],wt)
ggplot(FLQuants("WG"=m,"Lorenzen"=m.))+
  geom_line(aes(age,data,col=qname))

##Atlantic Skipjack
#Fonteneau & Pallarés (1999) assumed a constant M of 0.8,
skjPar=FLPar(c(linf=97.258,  k  =0.251,  t0=-0.01,
               a   =7.48e-6, b  =3.253,
               l50 =42,
               m1  =0.2795,    m2 =-0.288))
skjPar["a50"]=vonB(skjPar,length=skjPar["l50"])

m =FLQuant(rep(0.8,13),dimnames=list(age=0:12))
ln=vonB(skjPar,ages(m))
wt=len2wt(skjPar,ln)
m.=lorenzen(skjPar[c("m1","m2")],wt)

skjPar["m1"]=optimise(fn, c(0.001,10),m=m,par=skjPar[c("m1","m2")])$minimum

pars=FLPars(Bigeye  =betPar,
            Skipjack=skjPar,
            Yellowfin=yftPar)

save(pars,file="/home/laurie/Desktop/flr/git/FLife/data/pars.RData",compress="xz")

bet=gislasim(pars[["Bigeye"]])
bet=lh(bet,fnM=lorenzen)
ggplot(FLQuants("WG"=m,"Lorenzen"=m(bet)))+
  geom_line(aes(age,data,col=qname))

skj=gislasim(pars[["Skipjack"]])
skj=lh(skj,fnM=lorenzen)

yft=gislasim(pars[["Yellowfin"]])
yft=lh(yft,growth=gascuel,fnM=lorenzen)

eql=FLBRPs(Bigeye   =bet,
           Skipjack =skj,
           Yellowfin=yft)

save(eql,file="/home/laurie/Desktop/flr/git/FLife/data/eql.RData",compress="xz")
