## ----knitr_init, echo=FALSE, results="asis"------------------------------
library(knitr)

# output: rmarkdown::tufte_handout

## Global options
opts_chunk$set(echo    =TRUE,
               eval    =TRUE,
               cache   =!FALSE,
               cache.path="cache/",
               prompt  =FALSE,
               comment =NA,
               message =FALSE,
               tidy    =FALSE,
               warnings=FALSE,
               fig.height=4.5,
               fig.width =4.5)

## ----echo=TRUE,eval=TRUE,message=FALSE,warning=FALSE---------------------
library(ggplot2)
library(FLCore)
library(FLife)
library(ggplotFL)

## ----age-par-------------------------------------------------------------
par=lhPar(FLPar("linf"=100))
par

## ----age-----------------------------------------------------------------
age=FLQuant(0:20,dimnames=list(age=0:20))

ln =vonB(age,par)
wt =par["a"]*ln^par["b"]
wt =len2wt(ln,par)

## ----mat-----------------------------------------------------------------
mat=sigmoid(age,par)

## ---- fig.margin=TRUE, fig.cap="Biological properties",echo=FALSE--------
ggplot(FLQuants(Length=ln,Mass=wt,Maturity=mat))+
  geom_line(aes(age,data))+
  facet_grid(qname~.,scale="free")+
  xlab("Age")+ylab("")+theme_bw()

## ----age-m---------------------------------------------------------------
m =FLQuants(Lorenzen=lorenzen(wt,par[c("m1","m2")]),
            Gislason=gislason(ln,par[-(8:9)]))

## ---- fig.margin=TRUE, fig.height=4, fig.cap="Natural Mortality",echo=FALSE,eval=FALSE----
## ggplot(m)+
##   geom_line(aes(age,data,col=qname))+
##   xlab("Age")+ylab("Natural Mortality")+theme_bw()+
##   theme(legend.title=element_blank(),legend.position="bottom")

## ----age-sel-------------------------------------------------------------
sel=FLQuants("Flat"=dnormal(age,FLPar(a1="4",sl=3,sr=500)),
             "Dome"=dnormal(age,FLPar(a1="4",sl=3,sr=5)))

## ---- fig.margin=TRUE, fig.height=4, fig.cap="Selection patterns",echo=FALSE,eval=FALSE----
## ggplot(sel)+
##   geom_line(aes(age,data,col=qname))+
##   xlab("Age")+ylab("")+theme_bw()

## ----eval=FALSE,echo=FALSE-----------------------------------------------
## dat=data.frame(expand.grid(
## hg
