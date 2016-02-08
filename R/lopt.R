globalVariables(c("aaply"))


#' lopt
#'
#' @description Finds length at maximum biomass
#' 
#' @param params FLPar
#' @param ... any other arguments
#' 
#' @aliases lopt-method lopt,FLPar-method
#' 
#' @return FLPar with length at maximum biomass 
#' 
#' @details There are several ways to calculate \deqn{L_{opt}}, i.e.
#' i) \deqn{{2/3}^{rds}  L_{\infty}}
#' ii) \deqn{L_{\infty}\frac{3}{3+k/m}}
#' iii) by maximising the biomass of
#' iv) from an FLBRP object by fishing at F=0 and finding age where biomass is a maximum
#' 
#' @export
#' @docType methods
#' @rdname lopt
#' 
#' @seealso \code{\link{vonB}}  
#' 
#' @examples
#' \dontrun{
#' data(pars)
#' lopt(pars[[1]])
#' }
setGeneric('lopt', function(params,...)
  standardGeneric('lopt'))

setMethod("lopt", signature(params="FLPar"),
       function(params,
                   m=function(length,params) 
                       0.55*(length^-1.66)%*%(params["linf"]^1.44)%*%params["k"],
                   #m=function(length,params) params["m1"]%*%(exp(log(length)%*%params["m2"])),
                   ...){   
            
            dmns=dimnames(params)
            dmns$params="lopt"
            dm  =dim(params)
            
            res=aaply(params,seq(length(dm))[-1],function(x){
                   x.=FLPar(x)
                   rtn=try(optimise(loptFn,c(.01,c(x["linf"])*.99),params=x.,maximum=TRUE,m=m)$maximum)
                   if ("character" %in% mode(rtn)) rtn=NA
                   rtn})
            
            FLPar(array(res,dim=c(1, dm[-1]),dimnames=dmns))})

lopt_=function(params,
              m=function(length,params) 
                0.55*(length^-1.66)%*%(params["linf"]^1.44)%*%params["k"])
  
  optimise(loptFn,c(.01,c(params["linf"])*.99),params=params,maximum=TRUE,m=m)$maximum

loptFn=function(x,params,m,age=0:200){
  
  length=vonB(FLQuant(age,         dimnames=list(age=age)),params)
  m     =FLQuant(m(length,params),dimnames=list(age=age))
  mCum  =FLQuant(aaply(m,2,cumsum),dimnames=dimnames(m))
  
  a =qmax(vonB(params=params,length=FLQuant(x)),min(age))
  a =qmin(a,                                  max(age))
  
  a_=floor(a)
  
  if (is.na(a_)) return (as.numeric(NA))
  
  n =exp(-mCum[ac(a_)]-m[ac(a_)]*(a-a_))
  c(n*len2wt(x,params))}

#params=par[,10]
#params["t0"]=-.1
#loptFn(x,params)

loptFn2<-function(params) params["linf"]*3/(3+exp(params["m2"])/params["k"])

lopt3=function(x){
  
  fbar(x)          =fbar(x)[,1]*0
  stock.n(x)%*%stock.wt(x)}

#setMethod("lopt", signature(params="FLPar"),
    loptAge=function(params,
                   m     =function(length,params) params["m1"]%*%(exp(log(length)%*%params["m2"])),
                   growth=vonB,
                   ...){   

      loptFn=function(x,params,m){
        
        age   =0:ceiling(x)
        dmns  =list(age=age)
        length=vonB(age=FLQuant(pmin(age+0.5,x),dimnames=dmns),params=params)
        m.    =FLQuant(m(length,params),    dimnames=dmns)
        mCum  =FLQuant(aaply(m.,6,sum))
        n     =exp(-mCum)
        c(n*len2wt(length[ac(ceiling(x))],params))}
            
      dmns=dimnames(params)
      dmns$params="lopt"
      dm  =dim(params)
            
      res=aaply(params,seq(length(dm))[-1],function(x){
            x.=FLPar(x)
            rtn=try(optimise(loptFn,c(0,40),params=x.,maximum=TRUE,m=m)$maximum)
            if ("character" %in% mode(rtn)) rtn=NA
            rtn})
            
      vonB(FLQuant(FLPar(array(res,dim=c(1, dm[-1]),dimnames=dmns))),params)
      
      } #)

  