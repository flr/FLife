globalVariables(c("spr2v","srr2s"))


#' lhEql
#' 
#' Takes an \code{FLPar} object with life history and selectivity parameters
#' and generates an corresponding \code{FLBRP} object. Can uses a range of functional forms
#' 
#' @param params an \code{FLPar} object with life history parameters
#' @param ... any other arguments that include,
#' growth function for growth,
#' m function for natutal mortality,      
#' mat function for proportion mature-at-age,
#' sel function for selectivity-at-age,
#' sr character for stock recruitment relationship,
#' range age range,
#' spwn proportion of year when spawning occurrs, i.e. level of natural mortality prior to spawning,
#' fish proportion of year when fishing happens,
#' units units for FLQuant slots
#'  
#' @return \code{FLBRP} 
#' 
#' @export
#' @docType methods
#' @rdname lhEql
#' 
#' @aliases lhEql-method lhEql,FLPar-method
#' 
#' @seealso \code{\link{vonB}} \code{\link{lorenzen}} \code{\link{sigmoid}}  
#' 
#' @examples
#' \dontrun{
#' par=lhSim(FLPar(linf=100))
#' }
setGeneric('lhEql', function(params,...) standardGeneric('lhEql'))
setMethod("lhEql", signature(params='FLPar'),
          function(params,
            growth       =vonB,
            #fnM          =lorenzen, #function(params,length) exp(10e-3*params["m1"]%-%((log(length)%*%params["m2"]))),
            #fnM          =function(params,length,T=290,a=FLPar(c(a=-2.1104327,b=-1.7023068,c=1.5067827,d=0.9664798,e=763.5074169),iter=dims(params)$iter))
            #                       exp(a[1]%+%(a[2]%*%log(length/100))%+%(a[3]%*%log(params["linf"]/100))%+%(a[4]%*%log(params["k"]))%+%a[5]/T),
            m         =function(params,length) params["m1"]%*%(exp(log(length)%*%params["m2"])), 
            #fnM          =function(params,length) 0.55*(length^-1.61)%*%(params["linf"]^1.44)%*%params["k"],
            mat        =function(age,params) {
              a50=FLQuant(ceiling(rep(c(params["a50"]),each=dim(age)[1])),
                          dimnames=dimnames(age))
              res=FLQuant(0.5,dimnames=dimnames(age))
              res[age> a50]=1
              res[age< a50]=0
              res},
            sel          =dnormal,
            sr           ="bevholt",
            range        =c(min=0,max=40,minfbar=1,maxfbar=40,plusgroup=40),
            spwn         =c(params["a50"]-floor(params["a50"])),
            fish         = 0.5, # proportion of year when fishing happens
            units=if("units" %in% names(attributes(params))) attributes(params)$units else NULL,
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
                            iter=dimnames(params)$iter))
  # Get the lengths through different times of the year
  stocklen   <- growth(age+m.spwn,params) # stocklen is length at spawning time
  catchlen   <- growth(age+fish,  params) # catchlen is length when fishing happens
  
  midyearlen <- growth(age+0.5,params) # midyear length used for natural mortality
  
  # Corresponding weights
  swt=exp(log(stocklen%*%params["a"]))%*%params["b"]
  cwt=exp(log(catchlen%*%params["a"]))%*%params["b"]
  if ("bg" %in% dimnames(params)$param)  
    swt=exp(log(stocklen%*%params["a"]))%*%params["bg"]
  warning("FLPar%*%FLQuant operator sets 1st dim name to quant regardless")

  if ("numeric" %in% is(m)) m.=FLQuant(m,dimnames=dimnames(age)) else{
    if ("length" %in% names(formals(m)))   
      m.   =m(length=midyearlen,params=params) # natural mortality is always based on mid year length
    else if ("age" %in% names(formals(m))){ 
      m.   =m(age=age+0.5,params=params) # natural mortality is always based on mid year length
    }else if ("wt" %in% names(formals(m)))
      m.   =m(wt=swt,params[c("m1","m2")]) 
  
  names(dimnames(m.))[1]="age"}

  mat. =mat(age + m.spwn,params) # maturity is biological therefore + m.spwn 
  sel. =sel(age + fish,  params) # selectivty is fishery  based therefore + fish
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
 
  if (sr=="shepherd" & !("c" %in% dimnames(params)[[1]])){

    dmns=dimnames(params)
    
    dmns$params=c(dmns$params,"c")
    
    par.=FLPar(NA,dimnames=dmns)
    par.[dimnames(params)$params]=params
    par.["c"]=1
    params=par.}
  
  if (dims(params)$iter>1) {
    warning("Scarab, iters dont work for SRR:sv/ab etc")
    warning("Should be no need to specify mode of FLPar element") 
    
    if (sr=="shepherd")
      params(res)=FLPar(c(a=as.numeric(NA),b=as.numeric(NA),c=as.numeric(NA)),iter=dims(params)$iter)
    else
      params(res)=FLPar(c(a=as.numeric(NA),b=as.numeric(NA)),iter=dims(params)$iter)
    
    for (i in seq(dims(params)$iter))
      if (sr=="shepherd")
        params(res)[,i][]=unlist(c(FLCore::ab(params[c("s","v","c"),i],sr,spr0=FLCore::iter(spr0(res),i))[c("a","b","c")]))
    else
      params(res)[,i][]=unlist(c(FLCore::ab(params[c("s","v"),i],sr,spr0=FLCore::iter(spr0(res),i))[c("a","b")]))

    warning("iter(params(res),i)=ab(params[c(s,v),i],sr,spr0=iter(spr0(res),i))[c(a,b)] assignment doesnt work")
    warning("iter(FLBRP,i) doesn't work")
  }else{
    
    if (sr=="shepherd"){
      params(res)=FLCore::ab(params[c("s","v","c")],sr,spr0=spr0(res))[c("a","b","c")]
      }else{ 
      params(res)=FLCore::ab(params[c("s","v")],sr,spr0=spr0(res))[c("a","b")]
      }
    }

refpts(res)=propagate(FLBRP::refpts(res)[c("virgin","msy","crash","f0.1","fmax")],dims(params)$iter)
  res=brp(res)

  if ("fbar" %in% names(args)) 
    fbar(res)<-args[["fbar"]] else 
      if (any((!is.nan(FLBRP::refpts(res)["crash","harvest"])))) 
        fbar(res)<-FLQuant(seq(0,1,length.out=101),quant="age")%*%FLBRP::refpts(res)["crash","harvest"]
  
  names(dimnames(fbar(res)))[1]="age"
  res=brp(res)
   
  if (!("units" %in% names(attributes(params))))  return(res)
  if (all(is.na(attributes(params)$units)))  return(res)
  
  try(res <- setUnits(res, params),silent=TRUE)
  
  return(res)})
# 
# setMethod('ab', signature(x='FLparams', model='character'),
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

