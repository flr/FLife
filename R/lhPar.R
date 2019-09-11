utils::globalVariables(c("dlply"))
utils::globalVariables(c("mlply"))
utils::globalVariables(c("cast"))

lhValid=data.frame(old=c("linf","k","t0","a","b",
                         "ato95","a50","asym","bg","l50",
                         "s","v",
                         "a1","sl","sel2",
                         "m1","m2","m3"),
                    new=c("linf","k","t0","a","b",
                          "ato95","a50","asym","bg","l50",
                          "s","v",
                          "sel1","sel2","sel3",
                          "m1",  "m2",  "m3"),stringsAsFactors=FALSE)


#' @title Generates life history parameters
#' 
#' @description 
#' Uses life history theory to derive parameters for biological relationships, i.e. or growth, 
#' maturity, natural mortality. Selectivity by default is set so age at peak selectivity is the 
#' same as age at 50/% mature (a50) As a minimum all `lhPar` requires is `linf` the asymptotic 
#' length of the von Bertalannfy growth equation. 
#'  
#' @param   params \code{FLPar} object with parameters for life history equations and selection pattern.
#' Need Linfinity to estimate other parameters, if any other parameters supplied in \code{code} then
#' these are not provided by the algorithm 
#' @param   t0 of von Bertalanffy. This is a default that isnt normally derived
#' from life history theory, as are the following args.
#' @param a coefficient of length weight relationship
#' @param b exponent of length weight relationship
#' @param ato95 age at which 95\% of fish are mature, offset to age at which 50\% are mature
#' @param s steepness of stock recruitment relationship
#' @param v virgin biomass
#' @param sel1 selectivity-at-age parameter for double normal, age at maximum selectivity by default set to same as age at 100% mature
#' @param sel2 selectivity-at-age parameter for double normal, standard deviation of lefthand limb of double normal, by default 5
#' @param sel3 selectivity-at-age parameter for double normal, standard deviation of righthand limb of double normal, by default 5000
#' @param sl obsolete now replaced by sel2 
#' @param sr obsolete now replaced by sel3
#' @param m1 m-at-age parameter by default for Gislason empirical relationship
#' @param m2 m-at-age parameter, by default for Gislason empirical relationship
#' @param m3 m-at-age parameter, by default for Gislason empirical relationship
#' 
#' @export
#' 
#' @seealso \code{\link{loptAge}}, \code{\link{lhRef}}, \code{\link{lhPar}}, \code{\link{lhEql}}
#' 
#' @import methods
#' @docType methods
#' @rdname lhPar
#' @return object of class \code{FLPar} with missing parameters calculated from life history theory 
#' @examples
#' \dontrun{
#' lhPar(FLPar(linf=200))
#' }
#' 
# \deqn{ f(x) = \left\{
# \begin{array}{ll}
# 0 & x < 0 \\
# 1 & x \ge 0
# \end{array}
# \right. }{ (non-Latex version) }
lhPar=function(params,
               a=0.0003,b=3,
               ato95=1,
               sl=1,sr=5000,
               sel1=NA,sel2=1,sel3=5000,
               m1=0.55,m2=-1.61,m3=1.44,
               s=0.9,v=1000,
               t0=function(params)-exp(-0.3922-0.2752*log(params["linf"])%-%(1.038*log(params["k"])))
               ){
 
  if("data.frame"%in%class(params)) params=mf2FLPar(params)
 
  fn <- function(params,t0,a,b,ato95,sel2,sel3,s,v){
    
    names(dimnames(params)) <- tolower(names(dimnames(params)))
    
    if (!("a"     %in% dimnames(params)$params)) params=addpar(params,"a",   a)
    if (!("b"     %in% dimnames(params)$params)) params=addpar(params,"b",   b)
    if (!("bg"    %in% dimnames(params)$params)) {
      params=rbind(params,params["b"])
      dimnames(params)[[1]][length(dimnames(params)[[1]])]="bg"}
    if (!("s"     %in% dimnames(params)$params)) params=addpar(params,"s",   s)
    if (!("v"     %in% dimnames(params)$params)) params=addpar(params,"v",   v)
    
    ## growth parameters
    if (!("k"     %in% dimnames(params)$params)) {
      kpar  =FLPar(array(3.15*params["linf"]^(-0.64), dim=c(1, dims(params)$iter),dimnames=list(params="k", iter=seq(dims(params)$iter))))
  
      params=rbind(params,kpar) # From Gislason et al 2008, all species combined
      }
    
    if (!("t0"    %in% dimnames(params)$params)) params=addpar(params,"t0", t0(params))
    # Natural mortality parameters from Model 2, Table 1 Gislason 2010
    params=rbind(params,FLPar(m1=m1,m2=m2,m3=m3))
    
    #  if (!all(c("m1","m2")%in%dimnames(params)$params)){
    #    
    #     params=rbind(params,FLPar(m1= 0.55*(params["linf"]^1.44)%*%params["k"], iter=dims(params)$iter),
    #                         FLPar(m2=-1.61                                    , iter=dims(params)$iter))
    #    
    #    #m=(length^-1.61)%*%(exp(0.55)*params["linf"]^1.44)%*%params["k"]
    #    
    #    params=addpar(params,"m1", exp(0.55)*(params["linf"]^1.44)%*%params["k"])
    #    params=addpar(params,"m2", -1.61)
    #    }
  
    if (!("ato95" %in% dimnames(params)$params)) params=addpar(params,"ato95",ato95)  #rbind(params,FLPar("ato95" =ato95, iter=dims(params)$iter))
    if (!("sel2"    %in% dimnames(params)$params)) params=addpar(params,"sel2",   sel2)     #rbind(params,FLPar("sel2"    =sel2,    iter=dims(params)$iter))
    if (!("sel3"    %in% dimnames(params)$params)) params=addpar(params,"sel3",   sel3)     #rbind(params,FLPar("sel3"    =sel3,    iter=dims(params)$iter))
   
    ## maturity parameters from http://www.fishbase.org/manual/FishbaseThe_MATURITY_Table.htm
    if (!("asym"    %in% dimnames(params)$params)) params=params=addpar(params,"asym",  1) #rbind(params,FLPar("asym"    =asym, iter=dims(params)$iter))
    
    if (!("a50" %in% dimnames(params)$params)){
      if (!("l50" %in% dimnames(params)$params)){
        l50=0.72*params["linf"]^0.93
        dimnames(l50)$params="l50"
        }else{
        l50=params["l50"]
        }

      a50=log(1-(l50%/%params["linf"]))%/%(-params["k"])%+%params["t0"]
      dimnames(a50)$params="a50"
      
      params=rbind(params,a50)
      params=rbind(params,l50)
    }
    
    ## selectivity guestimate
    sel1=params["a50"]%+%params["ato95"]
    
    dimnames(sel1)$params="sel1"
   
    params=rbind(params,sel1)
    
    attributes(params)$units=c("cm","kg","1000s")
    
    return(params[lhValid$new[lhValid$new%in%dimnames(params)$params]])}
  
  if (dims(params)$iter==1) {
     return(fn(params[apply(is.na(params), 1, sum) == 0,],t0,a,b,ato95,sel2,sel3,s,v))
  }
  else{
    df=subset(as.data.frame(params),!is.na(data))
     
    res1=dlply(df,.(iter),function(x) as(x[,c("data","params")],"FLPar")[,1])
    
    res2=mlply(data.frame(iter=seq(length(res1))),function(iter,v1,v2,v3,v4,v5,v6,v7,v8) 
      fn(res1[[iter]],v1,v2,v3,v4,v5,v6,v7,v8),v1=t0,v2=a,v3=b,v4=ato95,v5=sel2,v6=sel3,v7=s,v8=v)
    
    res3=mdply(data.frame(iter=seq(length(res1))),function(iter) cbind(iter=iter,as.data.frame(iter(res2[[iter]],1))[,-2]))
    
    res4=as(res3,"FLPar")
    
   res4[]=unlist(c(cast(res3,params~iter,value="data")[,-1]))
    
    return(res4)}
  }

mf2FLPar=function(x){
  
  if ("iter"%in%names(x)){
     iters=x[,seq(length(dimnames(x)[[2]]))[dimnames(x)[[2]]=="iter"][1]]
     x    =x[,seq(length(dimnames(x)[[2]]))[dimnames(x)[[2]]!="iter"]]
  }else iter=seq(dim(x)[1])
  
  dmns=dimnames(x)[2:1]
  names(dmns)=c("params","iter")
  dmns[[2]]=seq(dim(x)[1])
  x=t(as.matrix(x))
  
  FLPar(array(x,dim=dim(x),dimnames=dmns),units="")}

addpar<-function(params,name,val)
  rbind(params,FLPar(array(val, dim=c(1, dims(params)$iter),dimnames=list(params=name, iter=seq(dims(params)$iter)))))


setUnits=function(res, par){

    if (is.null(attributes(params)$units)) return(res)
    units=attributes(params)$units
    allUnits=list("params"=      "",          
               "refpts"=         "",            
               "fbar"=           "",        
               "fbar.obs"=       "",    
               "landings.obs"=   paste(units[2],units[3]),
               "discards.obs"=   paste(units[2],units[3]),
               "rec.obs"=        units[3],         
               "ssb.obs"=        paste(units[2],units[3]),
               "stock.obs"=      paste(units[2],units[3]),
               "profit.obs"=     "",     
 #              "revenue.obs"=    "",    
               "landings.sel"=   "",    
               "discards.sel"=   "", 
               "bycatch.harvest"="",        
               "stock.wt"=       units[2],     
               "landings.wt"=    units[2],     
               "discards.wt"=    units[2],      
               "bycatch.wt"=     units[2],               
               "m"=              "",             
               "mat"=            "proportion", 
               "harvest.spwn"=   "proportion",          
               "m.spwn"=         "proportion",    
               "availability"=   "proportion",           
               "price"=          "",           
               "vcost"=          "",           
               "fcost"=          "")            

  
 res@units[names(allUnits)]=allUnits
    
 return(res)}
