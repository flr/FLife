#' @exportMethod FLStock
#' 
#' @title Generates an FLStock from life history parameters
#' 
#' @description 
#' Uses life history theory to derive parameters for biological relationships, i.e. or growth, 
#' maturity, natural mortality. Selectivity by default is set so age at peak selectivity is the 
#' same as age at 50\% mature (a50) As a minimum all `lhPar` requires is `linf` the asymptotic 
#' length of the von Bertalannfy growth equation. 
#'  
#' @param params \code{FLPar} object with parameters for life history equations and selection pattern.
#' Need Linfinity to estimate other parameters, if any other parameters supplied in \code{code} then
#' these are not provided by the algorithm 
#' @param t0 of von Bertalanffy. This is a default that isnt normally derived from life history theory, as are the following args.
#' @param a coefficient of length weight relationship
#' @param b exponent of length weight relationship
#' @param ato95 age at which 95\% of fish are mature, offset to age at which 50\% are mature
#' @param s steepness of stock recruitment relationship
#' @param v virgin biomass
#' @param sel1 selectivity-at-age parameter for double normal, age at maximum selectivity by default set to same as age at 100\% mature
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
#' }

lhStk<-function(...,
                k      =function(params,a=3.15,b=-0.64) a*params["linf"]^b,
                t0     =function(params,a=-0.3922,b=-0.2752,c=-1.038) 
                          -exp(a-b*log(params$linf)%-%(c*log(params$k))),
                l50    =function(params,a=0.72,b=0.93) a*params["linf"]^b,
                gowth  =vonB,
                mat    =logistic,
                sel    =dnormal,
                sr     ="bevholt",
                m      =list(model="gislason", params=c(m1=0.55, m2=-1.61, m3=1.44)),
                fmult  =function(x) refpts(x)["msy","harvest"]%*%FLQuant(seq(0,2,length.out=100)),
                range  =c(min=0,max=40,minfbar=1,maxfbar=40,plusgroup=40),
                spwn   =0, #c(params["a50"]-floor(params["a50"])),
                fish   =0.5, # proportion of year when fishing happens
                midyear=0.5){

  params=lhPar(...,m=m,k=k,t0=t0,l50=l50)

  eql   =lhEql(params,growth=vonB,
                      m      =if (is.FLQuant(m)) m else m$model,
                      sr     =sr,
                      mat    =mat,
                      sel    =sel,
                      range  =range,
                      spwn   =spwn,
                      fish   =fish,
                      midyear=midyear)
  
  fbar(eql)=fmult(eql)
  res=as(eql,"FLStock") 
  res=fwd(res,fbar=fbar(eql)[,-1],sr=eql)
  
  list(stk=res,eql=eql)}
  