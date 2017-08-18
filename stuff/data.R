library(FLCore)
load("/home/laurie/Desktop/MEGAsync/papers/fec/data/teleost.RData")
t.=as(teleost[,2:7],"FLPar")
dimnames(t.)$iter=teleost[,1]
attributes(t.)[["order"  ]]=as.factor(ac(teleost[,"order"]))
attributes(t.)[["genus"  ]]=as.factor(ac(teleost[,"genus"]))
attributes(t.)[["family" ]]=as.factor(ac(teleost[,"family"]))
attributes(t.)[["class." ]]=as.factor(ac(teleost[,"class"]))
attributes(t.)[["habitat"]]=as.factor(ac(teleost[,"habit"]))
teleost=t.

save(teleost,file="/home/laurie/Desktop/flr/FLife/data/teleost.RData",compress="xz")