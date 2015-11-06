#' lopt
#'
#' Finds length at maximum biomass
#' 
#' @param param
#' @param m natural mortality 
#' @param eql 
#' 
#' @return FLPar with length at maximum biomass 
#' 
#' @details There are several ways to calculate $L_{opt}$, i.e.
#' i) ${2/3}^{rds}$ of $L_{\infty}$
#' ii) $L_{\infty}\frac{3}{3+k/m}$
#' iii) by maximising the biomass of
#' iv) from an FLBRP object by fishing at F=0 and finding age where biomass is a maximum
#' 
#' 
#' @export
#' @docType methods
#' @rdname lopt
#' 
#' @seealso \code{\code{\link{vonB}}}  
#' 
#' @examples
#' \dontrun{
#' data(pars)
#' lopt(pars[[1]])
#' }
setGeneric('lopt', function(param,...)
  standardGeneric('lopt'))

setMethod("lopt", signature(param="FLPar"),
          function(param,
                   fnM=function(length,param) 
                     0.55*(length^-1.66)%*%(param["linf"]^1.44)%*%param["k"],...){   
            
            dmns=dimnames(param)
            dmns$params="lopt"
            dm  =dim(param)
            
            res=aaply(param,seq(length(dm))[-1],function(x){
                   x.=FLPar(x)
                   optimise(loptFn,c(.01,c(x["linf"])*.99),param=x.,maximum=TRUE,fnM=fnM)$maximum})
            
            FLPar(array(res,dim=c(1, dm[-1]),dimnames=dmns))})

lopt_=function(param,
              fnM=function(length,param) 
                0.55*(length^-1.66)%*%(param["linf"]^1.44)%*%param["k"])
  
  optimise(loptFn,c(.01,c(param["linf"])*.99),param=param,maximum=TRUE,fnM=fnM)$maximum


loptFn=function(x,param,fnM,age=0:200){
  
  length=vonB(FLQuant(age,         dimnames=list(age=age)),param)
  m     =FLQuant(fnM(length,param),dimnames=list(age=age))
  mCum  =FLQuant(aaply(m,2,cumsum),dimnames=dimnames(m))
  
  a =qmax(vonB(param=param,length=FLQuant(x)),min(age))
  a =qmin(a,                                  max(age))
  
  a_=floor(a)
  
  if (is.na(a_)) return (as.numeric(NA))
  
  n =exp(-mCum[ac(a_)]-m[ac(a_)]*(a-a_))
  c(n*len2wt(x,param))}

#params=par[,10]
#params["t0"]=-.1
#loptFn(x,params)

loptFn2<-function(params) params["linf"]*3/(3+exp(params["m2"])/params["k"])

lopt3=function(x){
  
  fbar(x)          =fbar(x)[,1]*0
  stock.n(x)%*%stock.wt(x)}

  