utils::globalVariables(c("rbind.fill","data","FLBRP2biodyn","quantity","optimise",
  "data","par","coefficients","lm","optimize","filter",
  "polygon","fft","lm","data","FLBRP2biodyn",
  "coefficients","data","fft","filter","polygon","quantity","rbind.fill",
  "dmns","rho","stk","aaply","maply","alply","laply","error",
  "%^%","%*%","%-%","%+%","pow",
  "FLPar",
  "lh","fbar<-","FLBRP","brp","refpts","refpts<-","as",
  "computeRefpts","computeRefpts","catch.obs",
  "fwdWindow","aaply","invGascuelFn","laply","lambda",
  "asym","ddply",".","year","a","b","grad","ddply",".","x",'aaply',"invgompertzFn",
  "spr2v","srr2s","fwd","jacobian","melt",
  "fwd","jacobian","melt","maply","jacobian","br","melt","maply","lambda","mdply","jacobian"))


setGeneric('knife', function(age,params,...) standardGeneric('knife'))
setGeneric('sigmoid', function(age,params,...) standardGeneric('sigmoid'))
setGeneric('dnormal', function(age,params,...) standardGeneric('dnormal'))
setGeneric('logistic', function(age,params,...) standardGeneric('logistic'))

setGeneric('gislason', function(length,params,...) standardGeneric('gislason'))
setGeneric('lorenzen', function(wt,params,...) standardGeneric('lorenzen'))

setGeneric('vonB',     function(age,params,...) standardGeneric('vonB'))
setGeneric('gascuel',  function(age,params,...) standardGeneric('gascuel'))
setGeneric('gompertz', function(age,params,...) standardGeneric('gompertz'))
setGeneric('richards', function(age,params,...) standardGeneric('richards'))

setGeneric('ages', function(object, ...) standardGeneric('ages'))
setGeneric('len2wt', function(length,params,...) standardGeneric('len2wt'))
setGeneric('wt2len', function(wt,params,...) standardGeneric('wt2len'))

setGeneric('mdd',   function(object,params,...) standardGeneric('mdd'))
setGeneric('matdd', function(age,params,...)    standardGeneric('matdd'))
setGeneric('grwdd', function(age,params,...)    standardGeneric('grwdd'))

setGeneric('lopt', function(params,...) standardGeneric('lopt'))
setGeneric('loptAge', function(params,...) standardGeneric('loptAge'))

setGeneric('leslie', function(object, ...) standardGeneric('leslie'))
#setGeneric('r', function(m,fec,...) standardGeneric('r'))
setGeneric('rod', function(object, ...) standardGeneric('rod'))

setGeneric('powh', function(len,n,...) standardGeneric('powh'))

setGeneric('rnoise',  function(n,len,...) standardGeneric('rnoise'))
setGeneric('rlnoise', function(n,len,...) standardGeneric('rlnoise'))

setGeneric('moment', function(object,...) standardGeneric('moment'))
setGeneric('sv', function(x,model, ...) standardGeneric('sv'))

setGeneric('lhEql', function(params,...) standardGeneric('lhEql'))

setGeneric('cc', function(age,n,...) standardGeneric('cc'))
