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


lhPar <- function(...,
    m=list(model="gislason", params=c(m1=0.55, m2=-1.61, m3=1.44))) {

  # DEFAULT defined parameter values + m
  args <- list(a=0.0003, b=3, ato95=1, sel2=1, sel3=5000,
    s=0.9, v=1000, asym=1, m=m)
  
  # PARSE ...
  input <- list(...)
  
  # FIND FLPar(s) and data.frame(s) in input
  flp <- unlist(lapply(input, function(x) is(x, "FLPar")))
  dtf <- unlist(lapply(input, function(x) is(x, "data.frame")))

  # CONVERT all to single list of vectors
  input <- c(input[!flp & !dtf],
      Reduce("c", lapply(input[flp], as, "list")),
      unlist(lapply(input[dtf], unlist)))

  # DROP any NAs
  input <- input[!unlist(lapply(input, function(x) any(is.na(x))))]

  # MERGE input and default args
  args[names(input)] <- input

  # PARSE m params
  margs <- args$m[unlist(lapply(args$m, is, "numeric"))][[1]]

  # ENSURE m params are named
  if(is.null(names(margs)))
    names(margs) <- paste0("m", seq(length(margs)))

  # STORE m model name
  mmodel <- args$m[unlist(lapply(args$m, is, "character"))][[1]]

  # ADD args and m params to form params
  params <- c(args[names(args) != "m"], margs)

  # EXPAND to max iters
  its <- max(unlist(lapply(params, length)))

  # CREATE output FLPar
  params <- do.call("FLPar", lapply(params, rep, length.out=its))
  
  # DERIVED parameters
  #
  # Gislason, H., J.G. Pope, J.C. Rice, and N. Daan. 2008. Coexistence in
  # North Sea fish communities: Implications for growth and natural mortality.
  # ICES J. Mar. Sci. 65 (4): 514–30.

  # k
  if(!"k" %in% dimnames(params)$params)
    params <- rbind(params, FLPar(k=3.15 * params$linf ^ (-0.64)))

  # t0
  if(!"t0" %in% dimnames(params)$params)
    params <- rbind(params, FLPar(t0=-exp(-0.3922 - 0.2752 *
      log(params$linf) %-% (1.038 * log(params$k)))))

  # l50 - a50
  if(!"l50" %in% dimnames(params)$params) {
    if("a50" %in% dimnames(params)$params) {
      params <- rbind(params,
        FLPar(l50=vonB(age=c(params$a50), params[c("k", "t0", "linf"),])))
    } else {
      params <- rbind(params, FLPar(l50=0.72 * params$linf ^ 0.93))
    }
  }
  if(!"a50" %in% dimnames(params)$params) {
    params <- rbind(params, FLPar(a50=log(1-(params$l50 %/%
      params$linf)) %/% (-params$k) %+% params$t0))
  }

  # sel1
  if(!"sel1" %in% dimnames(params)$params)
    params <- rbind(params, FLPar(sel1=params$a50 + params$ato95))
  
  # bg
  if(!"bg" %in% dimnames(params)$params)
    params <- rbind(params, FLPar(bg=params$b))

  # SORT params
  order <- c("linf", "l50", "a50", "ato95", "k", "t0", "a", "b",
    "m1", "m2", "m3", "sel1", "sel2", "sel3", "asym", "s", "v", "bg")

  params <- params[order,]

  # KEEP mmodel as attribute
  attr(params, "mmodel") <- mmodel

  return(params)
}


# ---

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
