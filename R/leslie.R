globalVariables(c("laply","lambda"))


#' leslie
#'
#' Creates a Leslie Matrix
#'  
#' @param object \code{FLBRP}
#' @param ... any other arguments
#' 
#' @aliases leslie-method leslie,FLBRP-method
#' 
#' @return \code{matrix} 
#' 
#' #' @export
#' @docType methods
#' @rdname leslie
#' 
#' @seealso \code{\link{lh}}  
#' 
#' @import FLBRP
#' 
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(FLCore)
#' library(FLBRP)
#' library(FLife)
#' data(pars)
#' pms=lhSim(pars[[1]])
#' eql=lh(pms,range = c(min=0,max=8, minfbar=1,maxfbar=8,plusgroup=8))
#' lsl=leslie(eql,fbar=c(refpts(eql)["crash","harvest"]))
#' }
setGeneric('leslie', function(object, ...)
  	standardGeneric('leslie'))

setMethod("leslie", signature(object="FLBRP"),
  function(object,fbar=FLQuant(0),numbers=TRUE,...){

  fbar(object)=fbar
  object=brp(object)
  names(dimnames(fbar(object)))[1]=names(dimnames(stock(object)))[1]
 
  ages=dims(object)$min:dims(object)$max
  
  mx=array(0, dim     =c(length(ages),length(ages),dims(stock(object))$iter),
              dimnames=list(age =ages,age=ages,
                            iter=seq(dims(stock(object))$iter)))
  
  #survivorship
  z=exp(-(m(object)))
  for (i in seq(dims(object)$iter)){
    diag(mx[-1,-length(ages),i]) =FLCore::iter(z[-length(ages)],i)
    if (range(object)["plusgroup"]==range(object)["max"])
      mx[length(ages),length(ages),i]=FLCore::iter(z[length(ages)],i)
    } 
  
  #recruitment
  tmp    = mat(object)*stock.wt(object)*stock.n(object)[,1]
  tmp2   = apply(tmp,2:6,sum)
  mx[1,,]= (rec(object)[,1]%*%tmp%/%tmp2)%/%stock.n(object)[,1]
  
  #Mass
  if (!numbers){
    #recruitment
    mx[1,,]=(rec(object)[,1]%*%stock.wt(object)[1,1]%/%ssb(object)[,1])%*%mat(object)[,1]
    
    #Growth
    incr=stock.wt(object)[-1,1]%/%stock.wt(object)[-length(ages),1]
    for (i in seq(dims(object)$iter))
      diag(mx[-1,-length(ages),i])=iter(incr[-length(ages)],i)*diag(mx[-1,-length(ages),i])              
    }

  #  diag(mx[-1,-length(ages)])=c(stock.wt(object)[-1,1])/c(stock.wt(object)[-length(ages),1])
  #  mx[1,]=c(stock.wt(object)[,1])*mx[1,]
  
  mx=FLPar(mx)
  mx[is.na(mx)]=0
  
  return(mx)})

#' r
#'
#' Estimates population growth rate for a Leslie matrix
#'  
#' @param m \code{FLPar}
#' @param fec \code{missing}
#' @param ... any other arguments
#' 
#' @aliases r-method r,FLPar-method
#' 
#' @return \code{FLPar} 
#' 
#' #' @export
#' @docType methods
#' @rdname lambda
#' 
#' @seealso \code{\link{leslie}}
#' 
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(FLBRP)
#' library(FLife)
#' data(pars)
#' pms=lhSim(pars[[1]])
#' eql=lh(pms)
#' lsl=leslie(eql)
#' lambda(lsl)
#' }
#setGeneric('r', function(m, ...) standardGeneric('r'))

setMethod("r", signature(m="FLPar",fec="missing"),
          function(m,...){

    object=m        
    #require(popbio)
    
    dmns=dimnames(object)[-2]
    dmns[1]="r"
    dm  =seq(length(dim(object)))[-(1:2)]
                        
    res=alply(object,dm,function(x) {
         rtn=try(lambda(x))
         if ("character" %in% mode(rtn)) rtn=NA
         rtn})
                        
    log(FLPar(array(res,dim     =unlist(laply(dmns,length)),
                        dimnames=dmns)))
    })
            
#setMethod("leslie", signature(object="FLBRP"),
oldLeslie=function(object,fbar=FLQuant(0),numbers=TRUE,...){
            
    args=list(...)  
    for (slt in names(args)[names(args) %in% names(getSlots("FLBRP"))])
        slot(object, slt)=args[[slt]]
    
    fbar(object)=fbar          
    ages=dims(object)$min:dims(object)$max
    mx  =matrix(0,nrow=length(ages),ncol=length(ages),dimnames=list(age=ages,age=ages))
            
    #survivorship
    diag(mx[-1,-length(ages)])   =exp(-m(object)[-length(ages)])
    if (range(object)["plusgroup"]==range(object)["max"])
      mx[length(ages),length(ages)]=exp(-m(object)[length(ages)])
            
    #recruitment
    tmp   = mat(object)*stock.wt(object)*stock.n(object)[,1]
    tmp2  = apply(tmp,2:6,sum)
    mx[1,]= rec(object)[,1]%*%tmp%/%tmp2%/%stock.n(object)[,1]
            
    if (!numbers){
      diag(mx[-1,-length(ages)])=diag(mx[-1,-length(ages)])*c(stock.wt(object)[-1,1])/c(stock.wt(object)[-length(ages),1])
      mx[1,]=c(stock.wt(object)[,1])*mx[1,]
      }
            
    mx[is.na(mx)]=0
            
    return(mx)}
#)



