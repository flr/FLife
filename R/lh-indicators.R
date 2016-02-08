globalVariables(c("a","b"))


maxAge=function(object,ref="msy",density=1/10000){
  
  fbar(object)=FLQuant(refpts(object)[ref,"harvest",drop=T])
  
  flag=stock.n(object)[,1]%/%stock.n(object)[1,1]<=density
  res=FLPar(array(0,c(1,dims(m(object))$iter),dimnames=list(param="maxage",iter=dimnames(m(object))$iter)))
  
  for (i in seq(dims(flag)$iter))
    iter(res,i)=min(ages(iter(flag,i))[iter(flag,i)])
  
  return(res)}

maxWt=function(object,ref="msy",density=1/10000)
  
   stock.wt[ac(maxAge(object,ref,density))]

maxLen=function(object,par,ref="msy",density=1/10000)
  
  (stock.wt[ac(maxAge(object,ref,density))]/a)^(1/b)
