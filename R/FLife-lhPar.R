# lhpar.R - DESC
# /lhpar.R
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the GPL 3.0

#' @title Generates life history parameters
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
#' #COMPARE with output of FLife::lhPar
#' 
#' x <- as(lhpar(linf=100), 'list')
#' x <- x[sort(names(x))]
#' y <- as(lhPar(FLPar(linf=100)), 'list')
#' y <- y[sort(names(y))]
#' 
#' all.equal(x,y)
#' 
#' for(i in seq(length(x)))
#'    cat(names(x[i]), ":", unlist(x[i]), "-", names(y[i]), ":", unlist(y[i]), "\n")
#'    
#'  # CALL with iters
#'  lhpar(FLPar(linf=100), v=rnorm(100, 300, 200))
#'  
#'  lhPar(FLPar(linf=rnorm(100, 80, 10)))
#'  lhPar(FLPar(linf=100, v=rnorm(100, 300, 200)))
#'  lhPar(FLPar(linf=100), FLPar(v=rnorm(100, 300, 200)))
#'  lhPar(FLPar(linf=100, v=rnorm(100, 300, 200)), t0=-1, data.frame(a=1,b=7))
#'  
#'  attributes(lhpar(FLPar(linf=100), v=rnorm(100, 300, 200)))$mmodel
#' }
#' 
# \deqn{ f(x) = \left\{
# \begin{array}{ll}
# 0 & x < 0 \\
# 1 & x \ge 0
# \end{array}
# \right. }{ (non-Latex version) }

lhPar <- function(...,
    m=list(model="gislason", params=c(m1=0.55, m2=-1.61, m3=1.44)),
    k  =function(params,a=3.15,b=-0.64) a*params["linf"]^b,
    t0 =function(params,a=-0.3922,b=-0.2752,c=-1.038) 
            -exp(a-b*log(params$linf)%-%(c*log(params$k))),
    l50=function(params,a=0.72,b=0.93) a*params["linf"]^b){
  
  args <- list(a=0.0003, b=3, bg=3, ato95=1, sel2=1, sel3=5000,
               s=0.9, v=1000, m=m, asym=1)
  
  # PARSE ...
  input <- list(...)

  # FIND FLPar(s) and data.frame(s)
  flp <- unlist(lapply(input, function(x) is(x, "FLPar")))
  dtf <- unlist(lapply(input, function(x) is(x, "data.frame")))

  # CONVERT all to single list
  input <- c(input[!flp & !dtf],
      Reduce("c", lapply(input[flp], as, "list")),
      unlist(lapply(input[dtf], unlist)))

  args[names(input)] <- input
      
  # PARSE m
  margs <- args$m[unlist(lapply(args$m, is, "numeric"))][[1]]

  if(is.null(names(margs)))
    names(margs) <- paste0("m", seq(length(margs)))

  mmodel <- args$m[unlist(lapply(args$m, is, "character"))][[1]]

  # CREATE params list
  params <- c(args[names(args) != "m"], margs)

  # SET iters
  its <- max(unlist(lapply(params, length)))

  # OUTPUT FLPar
  params <- do.call("FLPar", lapply(params, rep, length.out=its))
  # k
  if(!"k" %in% dimnames(params)$params)
    params <- rbind(params, FLPar("k"=k(params)))

  if(any(is.na(params["k"])))
    params["k",is.na(params["k"])]=FLPar("k"=k(params[,is.na(params["k"])]))

  # t0
  if(!"t0" %in% dimnames(params)$params)
    params <- rbind(params, FLPar("t0"=FLPar(t0(params))))
  if(any(is.na(params["t0"])))
    params["t0",is.na(params["t0"])]=FLPar("t0"=FLPar(t0(params[,is.na(params["t0"])])))
  
  # l50
  ## if l50 is NA and a50 supplied estimate l50
  if(("l50" %in% dimnames(params)$params)&("a50" %in% dimnames(params)$params)){
    flag=is.na(params["l50"]&!is.na(params["a50"]))
      if (any(flag))
        params["l50",flag]=vonB(c(params["a50",flag]),params[,flag])
    }

  ## if l50 not supplied but a50 is, then estimes l50
  if((!"l50" %in% dimnames(params)$params)&("a50" %in% dimnames(params)$params))
    params <- rbind(params, FLPar(a50=vonB(c(a50),params)))
  
  ## if l50 still missing estimate l50
  if(!"l50" %in% dimnames(params)$params)
    params <- rbind(params, FLPar("l50"=l50(params)))
  if(any(is.na(params["l50"])))
    params["l50",is.na(params["l50"])]=l50(params["l50",is.na(params["l50"])])
  
  # a50
  if(!"a50" %in% dimnames(params)$params)
    params <- rbind(params, FLPar(a50=NA))
  params["a50"]=log(1-(params$l50 %/%params$linf)) %/% (-params$k) %+% params$t0

  # sel1
  if(!"sel1" %in% dimnames(params)$params)
    params <- rbind(params, FLPar(sel1=params$a50 + params$ato95))
  
  # bg
  if(!"bg" %in% dimnames(params)$params)
    params <- rbind(params, FLPar(bg=params$b))

  # KEEP mmodel as attribute
  attr(params, "mmodel") <- mmodel

  return(params[c("linf","k","t0","a","b","l50","a50","ato95","asym","bg",
                  "m1","m2","m3","s","v","sel1","sel2","sel3")])}

