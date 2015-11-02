#' lh
#' 
#' Uses life history theory to derive parameters for biological relationships,
#' i.e. growth, maturity, natural mortality from.
#' 
#'
#' @param   \code{par} \code{FLPar} object with parameters for life history equations and selection pattern.
#' Need L_infty to estimate other parameters, if any other parameters supplied in \code{code} then
#' these are not provided by the algorithm 
#' @param   \code{t0} of von Bertalanffy. This is a default that isnt normally derived
#' from life history theory, as are the following args.
#' @param   \code{a} coefficient of length weight relationship
#' @param   \code{b} exponent of length weight relationship
#' @param   \code{ato95} age at which 95\% of fish are mature, offset to age at which 50\% are mature
#' @param   \code{sl} selectivity-at-age parameter, standard deviation of lefthand limb of double normal
#' @param   \code{sr} stock recruitment relationship
#' @param   \code{s} steepness of stock recruitment relationship
#' @param   \code{v} virgin biomass
#' 
#' @export
#' @docType methods
#' @rdname lh
#' @return An \code{FLPar} object with parameters
#' @examples
#' \dontrun{
#' lh(FLPar(linf=200))
#' }
#' 
# \deqn{ f(x) = \left\{
# \begin{array}{ll}
# 0 & x < 0 \\
# 1 & x \ge 0
# \end{array}
# \right. }{ (non-Latex version) }
lh=function(par,
            growth       =vonB,
            #fnM          =lorenzen, #function(par,len) exp(10e-3*par["m1"]%-%((log(len)%*%par["m2"]))),
            #fnM          =function(par,len,T=290,a=FLPar(c(a=-2.1104327,b=-1.7023068,c=1.5067827,d=0.9664798,e=763.5074169),iter=dims(par)$iter))
            #                       exp(a[1]%+%(a[2]%*%log(len/100))%+%(a[3]%*%log(par["linf"]/100))%+%(a[4]%*%log(par["k"]))%+%a[5]/T),
            fnM         =function(par,len) par["m1"]%*%(exp(log(len)%*%par["m2"])), 
            #fnM          =function(par,len) 0.55*(len^-1.61)%*%(par["linf"]^1.44)%*%par["k"],
            fnMat        =function(params,data) {
              a50=FLQuant(ceiling(rep(c(params["a50"]),each=dim(data)[1])),
                          dimnames=dimnames(data))
              res=FLQuant(0.5,dimnames=dimnames(data))
              res[data> a50]=1
              res[data< a50]=0
              res},
            fnSel        =dnormal,
            sr           ="bevholt",
            range        =c(min=1,max=40,minfbar=1,maxfbar=40,plusgroup=40),
            spwn         = 0,
            fish         = 0.5, # proportion of year when fishing happens
            units=if("units" %in% names(attributes(par))) attributes(par)$units else NULL,
            ...){

  # Check that spwn and fish are [0, 1]
  if (spwn > 1 | spwn < 0 | fish > 1 | fish < 0)
    stop("spwn and fish must be in the range 0 to 1\n")
  
  args<-list(...)

  if (("m.spwn" %in% names(args)))
    m.spwn =args[["m.spwn"]]
  else
    m.spwn=FLQuant(spwn, dimnames=list(age=range["min"]:range["max"]))
  
  if (("harvest.spwn" %in% names(args)))
    harvest.spwn =args[["harvest.spwn"]]
  else
    harvest.spwn=FLQuant(spwn, dimnames=list(age=range["min"]:range["max"]))
  
  age=FLQuant(range["min"]:range["max"],
              dimnames=list(age =range["min"]:range["max"],
                            iter=dimnames(par)$iter))
  # Get the lengths through different times of the year
  stocklen   <- growth(par,age+m.spwn) # stocklen is length at spawning time
  catchlen   <- growth(par, age+fish) # catchlen is length when fishing happens
  
  midyearlen <- growth(par, age+0.5) # midyear length used for natural mortality
  
  # Corresponding weights
  swt=exp(log(stocklen%*%par["a"]))%*%par["b"]
  cwt=exp(log(catchlen%*%par["a"]))%*%par["b"]
  if ("bg" %in% dimnames(par)$param)  
    swt=exp(log(stocklen%*%par["a"]))%*%par["bg"]
  warning("FLPar%*%FLQuant operator sets 1st dim name to quant regardless")

  if ("numeric" %in% is(fnM)) m.=FLQuant(fnM,dimnames=dimnames(age)) else{
    if ("len" %in% names(formals(fnM)))   
      m.   =fnM(par=par,len=midyearlen) # natural mortality is always based on mid year length
    else if ("age" %in% names(formals(fnM))){ 
      m.   =fnM(age=age+0.5,par=par) # natural mortality is always based on mid year length
    }else if ("wt" %in% names(formals(fnM)))
      m.   =fnM(par[c("m1","m2")],swt) 
  
  names(dimnames(m.))[1]="age"}

  #age<<-age
#mspwn<<-m.spwn
#return()
  mat. =fnMat(par,age + m.spwn) # maturity is biological therefore + m.spwn 
  sel. =fnSel(par,age + fish) # selectivty is fishery  based therefore + fish

  ## create a FLBRP object to   calculate expected equilibrium values and ref pts
  dms=dimnames(m.)

  res=FLBRP(stock.wt       =swt,
            landings.wt    =cwt,
            discards.wt    =cwt,
            bycatch.wt     =cwt,
            m              =m.,
            mat            =FLQuant(mat., dimnames=dimnames(m.)),
            landings.sel   =FLQuant(sel., dimnames=dimnames(m.)),
            discards.sel   =FLQuant(0,    dimnames=dimnames(m.)),
            bycatch.harvest=FLQuant(0,    dimnames=dimnames(m.)),
            harvest.spwn   =FLQuant(harvest.spwn,    dimnames=dimnames(m.)),
            m.spwn         =FLQuant(m.spwn,    dimnames=dimnames(m.)),
            availability   =FLQuant(1,    dimnames=dimnames(m.)),
            range          =range)
  ## FApex
  #if (!("range" %in% names(args))) range(res,c("minfbar","maxfbar"))[]<-as.numeric(dimnames(landings.sel(res)[landings.sel(res)==max(landings.sel(res))][1])$age)

  ## replace any slot passed in as an arg
  for (slt in names(args)[names(args) %in% names(getSlots("FLBRP"))[names(getSlots("FLBRP"))!="fbar"]])
    slot(res, slt)<-args[[slt]]

  params(res)=propagate(params(res),dims(res)$iter)
  ## Stock recruitment relationship
  model(res) =do.call(sr,list())$model
  
  if (sr=="shepherd" & !("c" %in% names(par))){
    
    dmns=dimnames(par)
    
    dmns$params=c(dmns$params,"c")
    
    par.=FLPar(NA,dimnamels=dmns)
    par.[dimnames(par)$params]=par
    par.["c"]=1
    par=par.}

  if (dims(par)$iter>1) {
    warning("Scarab, iters dont work for SRR:sv/ab etc")
    warning("Should be no need to specify mode of FLPar element") 
    
    params(res)=FLPar(c(a=as.numeric(NA),b=as.numeric(NA)),iter=dims(par)$iter)
    for (i in seq(dims(par)$iter))
      if (sr=="shepherd")
        params(res)[,i][]=unlist(c(FLCore:::ab(par[c("s","v","c"),i],sr,spr0=FLCore:::iter(spr0(res),i))[c("a","b","c")]))
    else
      params(res)[,i][]=unlist(c(FLCore:::ab(par[c("s","v"),i],sr,spr0=FLCore:::iter(spr0(res),i))[c("a","b")]))

    warning("iter(params(res),i)=ab(par[c(s,v),i],sr,spr0=iter(spr0(res),i))[c(a,b)] assignment doesnt work")
    warning("iter(FLBRP,i) doesn't work")
  }else{
    if (sr=="shepherd")
      params(res)=FLCore:::ab(par[c("s","v","c")],sr,spr0=spr0(res))[c("a","b","c")]
    else{ 
      params(res)=FLCore:::ab(par[c("s","v")],sr,spr0=spr0(res))[c("a","b")]
      }
    }

  refpts(res)=propagate(FLBRP:::refpts(res)[c("virgin","msy","crash","f0.1","fmax")],dims(par)$iter)
  res=brp(res)

  if ("fbar" %in% names(args)) 
    fbar(res)<-args[["fbar"]] else 
      if (any((!is.nan(FLBRP:::refpts(res)["crash","harvest"])))) 
        fbar(res)<-FLQuant(seq(0,1,length.out=101),quant="age")%*%FLBRP:::refpts(res)["crash","harvest"]
  
  names(dimnames(fbar(res)))[1]="age"
  res=brp(res)
   
  if (!("units" %in% names(attributes(par))))  return(res)
  if (all(is.na(attributes(par)$units)))  return(res)
  
  try(res <- setUnits(res, par),silent=TRUE)
  
  return(res)}
# 
# setMethod('ab', signature(x='FLPar', model='character'),
#   function(x, model, spr0=NA){
#  
#    s=x["a"]
#    v=x["b"]
#    a=FLPar(a=1,dimnames=dimnames(s))  
#    b=FLPar(b=1,dimnames=dimnames(v)) 
#    
#    if ("spr0" %in% dimnames(x)$params)
#       spr0=x["spr0"]  else 
#       spr0=FLPar(spr0,dimnames=dimnames(a)) 
# 
#    c=FLPar(c=1,dimnames=dimnames(a))  
#    d=FLPar(d=1,dimnames=dimnames(a))  
#    if (("c" %in% dimnames(x)$params))  c=x["c"]
#    if (("d" %in% dimnames(x)$params))  d=x["d"]
# 
#    v <- v*spr2v(model, spr0, a, b, c, d)
#    s <- s*srr2s(model, ssb=v*.2, a=a, b=b, c=c, d=d) / srr2s(model, ssb=v, a=a, b=b, c=c, d=d)
#   
#    res=rbind(s, v, spr0)
#  
#    if ("c" %in% dimnames(x)$params)
#      res=rbind(res, c)
#  
#    if ("d" %in% dimnames(x)$params)
#      res=rbind(res, d)
#  
#    res=rbind(res, spr0)
#  
#    return(res)})

