stats<-function(x){
  
  mn=aaply(x,1,mean, na.rm=TRUE)
  sd=aaply(x,1,var,  na.rm=TRUE)^0.5
  n =aaply(x,1,function(x) sum(!is.na(x)))
  se=sd/n^0.5
  
  data.frame(param=names(se),cbind(mn=mn,sd=sd,se=se,n=n))}

sim<-function(x,niters=500,se=0.3){
  
  mn=aaply(x,1,mean, na.rm=TRUE)
  sd=aaply(x,1,var,  na.rm=TRUE)^0.5
  n =aaply(x,1,function(x) sum(!is.na(x)))
  se=sd/n^0.5
  
  if (any(is.na(se))) se[is.na(se)]=se
  
  y=data.frame(mn=mn,se=se)
  y=mdply(y,function(mn,se) rnorm(niters,mn,se))[,-(1:2)]
  
  res=FLPar(array(unlist(c(y)),c(dim(x)[1],niters)))
  
  dimnames(res)$params=names(mn)
  
  res}
