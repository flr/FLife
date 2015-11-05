#' lopt
#'
#' Finds length at maximum biomass
#' 
#' @param param
#' @param fnM natural mortality 
#' 
#' @return numeric with length 
#' 
#' #' @export
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
            
            optimise(loptFn,c(.01,c(param["linf"])*.99),param=param,maximum=TRUE,fnM=fnM)$maximum
            })

lopt_=function(param,
              fnM=function(length,param) 
                0.55*(length^-1.66)%*%(param["linf"]^1.44)%*%param["k"])
  
  optimise(loptFn,c(.01,c(param["linf"])*.99),param=param,maximum=TRUE,fnM=fnM)$maximum


loptFn=function(x,len,age=0:200,
                param,fnM){
  
  length=vonB(FLQuant(age,         dimnames=list(age=age)),param)
  m     =FLQuant(fnM(length,param),dimnames=list(age=age))
  mCum  =FLQuant(aaply(m,2,cumsum),dimnames=dimnames(m))
  
  a =max(0,vonB(param=param,length=FLQuant(x)))
  
  a_=floor(a)
  
  n =exp(-mCum[ac(a_)]-m[ac(a_)]*(a-a_))
  c(n*len2wt(x,param))}

#params=par[,10]
#params["t0"]=-.1
#loptFn(x,params)

loptFn2<-function(params) params["linf"]*3/(3+exp(params["m2"])/params["k"])

lopt3=function(x){
  
  fbar(x)          =fbar(x)[,1]*0
  stock.n(x)%*%stock.wt(x)}

  