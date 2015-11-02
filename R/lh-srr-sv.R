#### SRR 
setGeneric('sv', function(x,model, ...)
  standardGeneric('sv'))

setMethod('sv', signature(x='FLPar', model='character'),
  function(x, model, spr0=NA){
 
   a=x["a"]
   b=x["b"]
   s=FLPar(a=1,dimnames=dimnames(a))  
   v=FLPar(b=1,dimnames=dimnames(a))  
   spr0=FLPar(spr0,dimnames=dimnames(a))  

   if ("spr0" %in% dimnames(x)$params)
     spr0=x["spr0"] 

   c=FLPar(c=1,dimnames=dimnames(a))  
   d=FLPar(d=1,dimnames=dimnames(a))  
   if (("c" %in% dimnames(x)$params))  c=x["c"]
   if (("d" %in% dimnames(x)$params))  d=x["d"]

   v <- v*spr2v(model, spr0, a, b, c, d)
   s <- s*srr2s(model, ssb=v*.2, a=a, b=b, c=c, d=d) / srr2s(model, ssb=v, a=a, b=b, c=c, d=d)
  
   res=rbind(s, v, spr0)
 
   if ("c" %in% dimnames(x)$params)
     res=rbind(res, c)
 
   if ("d" %in% dimnames(x)$params)
     res=rbind(res, d)
 
   res=rbind(res, spr0)
 
   return(res)})

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
