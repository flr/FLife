knifeEdgeMsy=function(obj){
  
  fbar(obj)=FLQuant(0,dimnames=list(age=0,year=1))
  
  as.data.frame((stock.n(obj)*stock.wt(obj)),drop=TRUE)}