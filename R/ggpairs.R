utils::globalVariables(c("ggplot","geom_bar","ccf","geom_density","geom_point","geom_smooth"))

my_density<-function(data,mapping,...){
  ggplot(data=data,mapping=mapping)+
    geom_density(...,lwd=1)}

my_bar <- function(data,mapping,...){
  ggplot(data=data,mapping=mapping)+
    geom_bar(...)}

my_smooth <- function(data,mapping,...){
  ggplot(data=data,mapping=mapping)+
    geom_point(...,size=.5)+
    geom_smooth(...,method="lm",se=FALSE)}

my_ccf<-function(cpue){
  cc=mdply(expand.grid(a=names(cpue),b=names(cpue)),
           function(a,b){
             #print(paste(a,b))
             res=model.frame(mcf(FLQuants(cpue[c(a,b)])))
             res=subset(res,!is.na(res[,7])&!is.na(res[,8]))
             
             if (dim(res)[1]>10){
               res=data.frame(lag=-7:7,data=ccf(res[,7],res[,8],plot=F,
                                                lag.max=7)$acf)
               return(res)}else{return(NULL)}}
  )}

