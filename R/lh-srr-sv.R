spr2v <- function(model, spr, a=NULL, b=NULL, c=NULL, d=NULL){
  # SSB as function of ssb/rec
  return(switch(model,
                "bevholt"  = a*(spr)-b,
                "ricker"   = log(a*spr)/b,
                "cushing"  = (1/(a*spr))^(1/(b-1)),
                "shepherd" = b*(a*spr-1)^(1/c),
                "segreg"   = ifelse(ssb <= b, a/(spr), 0),
                "mean"     = a/(spr),
                "dersh"    = ssb*a*(1-b*c*ssb)^c,
                "pellat"   = 1/(a/ssb-a/ssb*(ssb/b)^c),
                NULL))}

srr2s <- function(model, ssb=NULL, spr=NULL, a=NULL, b=NULL, c=1, d=NULL)
{
  #recruits as function of ssb or ssb/rec
  if (is.null(ssb) & !is.null(spr))
    ssb <- spr2v(model, spr, a, b, c, d)
  
  eval(as.list(do.call(model, list())$model)[[3]], envir=list(ssb=ssb, spr0=spr, a=a, b=b, c=c, d=d))
} # }}}


#' @title Calculates steepness and virgin biomass
#' 
#' @description 
#' Calculates steepness and virgin biomass given a and b for a Beverton and Holt SRR
#' 
#' @param x \code{FLPar} with a and b
#' @param model \code{character} with name of stock recruitment relationship, by default "bevholt"
#' @param spr0 \code{} spawner per recruit at F=0
#' @param ... any other arguments
#'  
#' @return  \code{FLPar} with values for steepness (s) and virgin biomass (v) 
#' @export
#' @docType methods
#' @rdname sv
#' 
#' 
#' @examples
#' \dontrun{
#' #bug
#' params=FLPar(a=37.8,b=8.93)
#' sv(params,"bevholt",.4)
#' }
#
setMethod('sv', signature(x='FLPar', model='character'),
  function(x, model, spr0=NA){
 
   a=x["a"]
   b=x["b"]
   dmns=dimnames(x["a"])
   dmns[[1]]="s"
   s=FLPar(1,dimnames=dmns)
   dmns[[1]]="v"
   v=FLPar(1,dimnames=dmns)  
   dmns[[1]]="spr0"
   spr0=FLPar(spr0[drop=T],dimnames=dmns)  

   if ("spr0" %in% dimnames(x)$params)
     spr0=x["spr0"] 

   c=FLPar(1,dimnames=dimnames(a))  
   d=FLPar(1,dimnames=dimnames(a))  
   if (("c" %in% dimnames(x)$params))  c=x["c"]
   if (("d" %in% dimnames(x)$params))  d=x["d"]

   v <- v*spr2v(model, spr0, a, b, c, d)
   s <- s*srr2s(model, ssb=v*.2, a=a, b=b, c=c, d=d) / srr2s(model, ssb=v, a=a, b=b, c=c, d=d)
  
   res=rbind(s, v, spr0)
 
   if ("c" %in% dimnames(x)$params)
     res=rbind(res, c)
 
   if ("d" %in% dimnames(x)$params)
     res=rbind(res, d)
 
   return(res)})

#x=FLPar(s=0.75,v=1500,spr0=12)
#sv(ab(x,"bevholt"),"bevholt")

abPars. <- function(x,spr0=NA,model){
  s=x["s"]
  v=x["v"]
  if ("c" %in% names(x))
     c=x["c"]
  if ("d" %in% names(x))
     d=x["d"]
  if ("spr0" %in% names(x))
     spr0=x["spr0"]
  # converts a & b parameterisation into steepness & virgin biomass (s & v)
  switch(model,
    "bevholt"   ={a=(v%+%(v%-%s%*%v)%/%(5%*%s%-%1))%/%spr0; b=(v%-%s%*%v)%/%(5%*%s%-%1)},
    "bevholtSV" ={a=(v+(v-s*v)/(5*s-1))/spr0; b=(v-s*v)/(5*s-1)},
    "ricker"    ={b=log(5*s)/(v*0.8); a=exp(v*b)/spr0},
    "rickerSV"  ={b=log(5*s)/(v*0.8); a=exp(v*b)/spr0},
    "cushing"   ={b=log(s)/log(0.2); a=(v^(1-b))/(spr0)},
    "cushingSV" ={b=log(s)/log(0.2); a=(v^(1-b))/(spr0)},
    "shepherd"  ={b=v*(((0.2-s)/(s*0.2^c-0.2))^-(1/c)); a=((v/b)^c+1)/spr0},
    "shepherdSV"={b=v*(((0.2-s)/(s*0.2^c-0.2))^-(1/c)); a=((v/b)^c+1)/spr0},
    "mean"      ={a=v/spr0;b=NULL},
    "meanSV"    ={a=v/spr0;b=NULL},
    "segreg"    ={a=5*s/spr0; b=v/(a*spr0)},
    "segregSV"  ={a=5*s/spr0; b=v/(a*spr0)},
    {stop("model name not recognized")})

  res <- c(a=a, b=b)
  return(res[!is.null(res)])}
