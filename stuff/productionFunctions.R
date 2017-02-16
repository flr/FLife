library(FLBRP)
library(FLife)
  
data(teleost)

alb=lhPar(teleost[,"Thunnus alalunga"])
msy=lhEql(alb)

#ref=refpts(msy)[c("msy","virgin"),c("harvest","ssb","yield")]

par  =FLPar(m1=3.0,m2=-0.288)
refY =stock.n(msy)%*%stock.wt(msy)
scale=(stock.n(eql)[,1]%*%stock.wt(eql)%-%refY)%/%refY

virgin=lhEql(alb,m=c(mdd(stock.wt(virgin),par,scale)))

computeRefpts(msy)
computeRefpts(virgin)
