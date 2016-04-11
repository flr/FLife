fapexAge<-function(object){
  tmp=harvest(object)
  tmp[]=fapex(object)
  tmp=FLQuant(ages(tmp)[harvest(object)==tmp],
              dimnames=dimnames(fapex(object)))
  tmp}