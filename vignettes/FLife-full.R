## ----knitr_init, echo=FALSE, results="hide"------------------------------
library(knitr)
## Global options
opts_chunk$set(cache     =TRUE,
               echo      =!TRUE,
               eval      =TRUE,
               prompt    =FALSE,
               comment   =NA,
               message   =FALSE,
               warning   =FALSE,
               tidy      =FALSE,
               fig.height=6,
               fig.width =8)

iFig=0


## ------------------------------------------------------------------------
library(ggplot2)
library(reshape)
library(plyr)


## ----install,echo=TRUE,eval=FALSE----------------------------------------
## install.packages("FLife", repos = "http://flr-project.org/R")


## ----install-64,echo=FALSE,eval=FALSE------------------------------------
## install.packages("FLife", repos = "http://flr-project.org/R", INSTALL_opts="-no_multiarch" )


## ----lib,echo=TRUE-------------------------------------------------------
library(FLife)


## ----data,echo=TRUE------------------------------------------------------
data(teleost)


## ----data-3,echo=TRUE----------------------------------------------------
teleost


## ---- fig.height=8,fig.cap="Relationship between life history parameters."----
library(GGally)

habitat=ifelse(attributes(teleost)$habitat=="demersal","Demersal","Other")

my_smooth <- function(data,mapping,...){
  ggplot(data=data,mapping=mapping)+
  geom_point(...,size=.5)+
  geom_smooth(...,method="lm",se=FALSE)}

my_density <- function(data,mapping,...){
  ggplot(data=data,mapping=mapping)+
  geom_density(...,lwd=1)}
  
ggpairs(cbind(transform(model.frame(teleost)[,-c(7)],
                        linf=log(linf),k=log(k),l50=log(l50))),
  lower = list(continuous = wrap(my_smooth)),
  diag=list(continuous=wrap(my_density,alpha=0.2)),
  title = "")+
  theme(legend.position ="bottom",
  panel.grid.major =element_blank(),
  axis.ticks       =element_blank(),
  axis.text.x      =element_text(angle=-45),
  axis.text.y      =element_blank(),
  panel.border     =element_rect(linetype = 1, colour="black", fill=NA))+
  theme_bw()


## ----echo=TRUE,fig.cap="Relationship between k and Linfinty"-------------
linf=teleost["linf"]

hat=lhPar(linf)

ggplot()+
  geom_line(aes(linf,k), data=model.frame(hat))+
  geom_point(aes(linf,k),data=model.frame(teleost),col="red")+
  scale_x_log10()+
  scale_y_log10()+
  xlab(expression(L[infinty]))


## ------------------------------------------------------------------------
lhPar(linf)


## ---- fig.caption="Estimates of steepness"-------------------------------
y=2.706-3.698*teleost["l50"]/teleost["linf"]

invLogit<-function(y) 0.2+exp(y)/(1+exp(y))

Steepness=invLogit(y)

hist(Steepness)


## ---- devtools, echo=TRUE, eval=FALSE------------------------------------
## 	library(devtools)
## 	install_github("flr/FLife")

