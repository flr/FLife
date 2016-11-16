qnt <-function(x) rank(x)/length(x)

comp<-function(x) {
  res=apply(x,2,qnt)
  apply(res,1,mean)}

pow<-function(a,b) a^b
