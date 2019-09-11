utils::globalVariables(c("len"))

setGeneric('m', function(object,model,params,...) 
  standardGeneric('m'))
setGeneric('gislason', function(length,params,...) 
  standardGeneric('gislason'))
setGeneric('lorenzen', function(wt,params,...) 
  standardGeneric('lorenzen'))

mNms<-c("gislason",
        "lorenzen",
        "rikhter",
        "rikhter2",
        "jensen",
        "jensen2",
        "griffiths",
        "chen",
        "roff",
        "petersen",
        "djababli",
        "jensen2",
        "charnov")

mPar<-function(model){
  
  model=as.character(model)

  res=array(c(
    c(0.55,-1.61, 1.44),
    c(1.521,0.72,-0.155),
    c(1.65, NA,   NA),
    c(NA,   NA,   NA),
    c(1.5,  NA,   NA),
    c(NA,   NA,   NA),
    c(1.406,-0.096,0.78),
    c(1.0661,-0.1172,0.5092),
    c(3.0,    NA,    NA),
    c(1.28,-0.25,    NA),
    c(NA,     NA,    NA),
    c(NA,     NA,    NA),
    c(NA,     NA,    NA),
    c(NA,     NA,    NA)),c(3,13,1))
  
  res=aperm(res,c(2,1,3))
  
  res=FLPar(aperm(res,c(2,1,3)),
        dimnames=list(params=c("m1","m2","m3"),
                      model=mNms,
                      iter=1))
  
  FLPar(res[,model,drop=T])}

mVar       =c("age",  "age",   "age",     "age", "age",   "age",    "len",    "len",     "age",     "wt")
names(mVar)=c("roff","rikhter","rikhter2","chen","jensen","jensen2","charnov","gislason","petersen","lorenzen") 


mFn<-function(model,flq,params){
  
  switch(as.character(model),

  gislason={
      exp(params["m1"]%-%(params["m1"]%*%log(flq))+(params["m1"]%*%log(params["linf"]))%+%log(params["k"]))},
  
  roff={
    res=(3*params["k"]%*%params["linf"])*(1.0-params["l50"]%/%params["linf"])%/%params["l50"]
    
    dimnames(res)$params="m"
    res},
  rikhter={
    tm=params["t0"]%+%log(1-params["l50"]/(params["linf"]))%/%(-params["k"])
    res=params["b"]%*%params["k"]%/%(exp(params["k"]%*%(tm%-%params["t0"]))-1)
    
    dimnames(res)$params="m"
    res},
  
  rikhter2={
    tm=params["t0"]%+%log(1-params["l50"]/(params["linf"]))%/%(-params["k"])
    res=1.521/tm^0.73-0.155
    
    dimnames(res)$params="m"
    res},
  
  griffiths={
    winf=params["a"]%*%(params["linf"]^params["b"])
    res=(1.406*(winf^-0.096))%*%params["k"]^0.78
    
    dimnames(res)$params="m"
    res},
  
  djababli={
    res=(1.066*params["linf"]^-0.1172)%*%params["k"]^0.5092
    
    dimnames(res)$params="m"
    res},
  
  jensen={
    res=1.5*params["k"]
    
    dimnames(res)$params="m"
    res},
  
  jensen2={
    tm=params["t0"]%+%log(1-params["l50"]/(params["linf"]))%/%(-params["k"])
    res=1.65/tm
    
    dimnames(res)$params="m"
    res},
  
  charnov={
    res=params["k"]%*%(params["linf"]%/%len)^1.5
    
    res},
  
  petersen={
    1.28*wt^(-0.25)},
  
  chen={
    m =params["k"]/(1-exp(-params["k"]%*%(age%-%params["t0"])))
    
    tm =-(1/params["k"])*log(1-exp(params["k"]*params["t0"]))+params["t0"]
    bit=exp(-params["k"]*(tm-params["t0"]))
    
    a0=1-bit
    a1=params["k"]*bit
    a2=-0.5*params["k"]^2*bit
    age.=age>c(tm)
    m[age.] =params["k"]/(a0+a1*(age[age.]-tm)+a2*(age[age.]-tm)^2)
    
    dimnames(m)$params="m"
    return(m)}) 
  }

setMethod('m', signature(object='FLQuant',model="character",params='FLPar'),
          function(object,model,params,...) { 

          mFn(model,object,params)})

setMethod('m', signature(object='FLStock',model="character",params='FLPar'),
          function(object,model,params,...) { 

          flq=switch(mVar[model],
                 len=wt2len(stock.wt(object),params),
                 wt =       stock.wt(object),
                 age=  ages(stock.wt(object))) 
                 
          m(flq,model,params)})

setMethod('m', signature(object='FLBRP',model="character",params='FLPar'),
          function(object,model,params,...) { 

            flq=switch(mVar[model],
                       len=wt2len(stock.wt(object),params),
                       wt =       stock.wt(object),
                       age=  ages(stock.wt(object))) 
            
            m(flq,model,params)})
setMethod('m', signature(object='FLQuant',model="factor",params='FLPar'),
          function(object,model,params,...) { 
          
          model=as.character(model)
            
          m(flq,model,params)})

setMethod('m', signature(object='FLStock',model="factor",params='FLPar'),
          function(object,model,params,...) { 
            
            model=as.character(model)
            
            flq=switch(mVar[model],
                       len=wt2len(stock.wt(object),params),
                       wt =       stock.wt(object),
                       age=  ages(stock.wt(object))) 
            
            m(flq,model,params)})
setMethod('m', signature(object='FLBRP',model="factor",params='FLPar'),
          function(object,model,params,...) { 
            model=as.character(model)
            flq=switch(mVar[model],
                       len=wt2len(stock.wt(object),params),
                       wt =       stock.wt(object),
                       age=  ages(stock.wt(object))) 
            
            m(flq,model,params)})
setMethod('m', signature(object='FLQuant',model="character",params='numeric'),
          function(object,model,params,...) { 
          
          params=FLPar(params)    
          m(flq,model,params)})
setMethod('m', signature(object='FLStock',model="character",params='numeric'),
          function(object,model,params,...) { 
            
            params=FLPar(params)    
            
            flq=switch(mVar[model],
                       len=wt2len(stock.wt(object),params),
                       wt =       stock.wt(object),
                       age=  ages(stock.wt(object))) 
            
            m(flq,model,params)})
setMethod('m', signature(object='FLBRP',model="character",params='numeric'),
          function(object,model,params,...) { 
            
            params=FLPar(params)    
            
            flq=switch(mVar[model],
                       len=wt2len(stock.wt(object),params),
                       wt =       stock.wt(object),
                       age=  ages(stock.wt(object))) 
            
            m(flq,model,params)})
setMethod('m', signature(object='FLQuant',model="factor",params='numeric'),
          function(object,model,params,...) { 
          
          model=as.character(model)
          params=FLPar(params)
          
          m(flq,model,params)})
setMethod('m', signature(object='FLStock',model="factor",params='numeric'),
          function(object,model,params,...) { 
          
            model=as.character(model)
            params=FLPar(params)    
            
            flq=switch(mVar[model],
                       len=wt2len(stock.wt(object),params),
                       wt =       stock.wt(object),
                       age=  ages(stock.wt(object))) 
            
            m(flq,model,params)})
setMethod('m', signature(object='FLBRP',model="factor",params='numeric'),
          function(object,model,params,...) { 
            
            model=as.character(model)
            params=FLPar(params)
            
            flq=switch(mVar[model],
                       len=wt2len(stock.wt(object),params),
                       wt =       stock.wt(object),
                       age=  ages(stock.wt(object))) 
            
            m(flq,model,params)})
setMethod('m', signature(object='FLQuant',model="character",params='missing'),
          function(object,model,params,...) { 
            
            params=mPar(model)
            
            m(flq,model,params)})
setMethod('m', signature(object='FLStock',model="character",params='missing'),
          function(object,model,params=mPar(model),...) { 
            
            params=mPar(model)
            
            flq=switch(mVar[model],
                       len=wt2len(stock.wt(object),params),
                       wt =       stock.wt(object),
                       age=  ages(stock.wt(object))) 
            
            m(flq,model,params)})
setMethod('m', signature(object='FLBRP',model="character",params='missing'),
          function(object,model,params,...) { 
            
            params=mPar(model)
            
            flq=switch(mVar[model],
                       len=wt2len(stock.wt(object),params),
                       wt =       stock.wt(object),
                       age=  ages(stock.wt(object))) 

            m(flq,model,params)})
setMethod('m', signature(object='FLQuant',model="factor",params='missing'),
          function(object,model,params,...) { 
          
          params=mPar(model)
            
          model=as.character(model)
          m(flq,model,params)})

setMethod('m', signature(object='FLStock',model="factor",params='missing'),
          function(object,model,params,...) { 
          
            params=mPar(model)
            
            model=as.character(model)
            
            flq=switch(mVar[model],
                       len=wt2len(stock.wt(object),params),
                       wt =       stock.wt(object),
                       age=  ages(stock.wt(object))) 
            
            m(flq,model,params)})
setMethod('m', signature(object='FLBRP',model="factor",params='missing'),
          function(object,model,params,...) { 
            
            params=mPar(model)
            
            model=as.character(model)
            params=mPar(model)
            
            flq=switch(mVar[model],
                       len=wt2len(stock.wt(object),params),
                       wt =       stock.wt(object),
                       age=  ages(stock.wt(object))) 
            
            m(flq,model,params)})

#' gislason
#' @description 
#' gislason natural mortality relatoinship estimate M as a function of length. 
#' 
#' @param length at which M is to be predicted
#' @param m1 0.55
#' @param m2 1.44
#' @param m3 -1.61
#' @param ... any other arguments
#' 
#' @aliases  gislason gislason-method 
#'           gislason,FLQuant,FLPar-method 
#'           gislason,FLQuant,missing-method 
#'           gislason,FLQuant,numeric-method
#'
#' @export
#' @docType methods
#' @rdname gislason
#' 
#' @seealso \code{\link{lorenzen}}
#'  
#' @examples
#' \dontrun{
#' params=lhPar(FLPar(linf=111))
#' len=FLQuant(c( 1.90, 4.23, 7.47,11.48,16.04,20.96,26.07,31.22,
#'                36.28,41.17,45.83,50.20,54.27,58.03,61.48,64.62),
#'              dimnames=list(age=1:16))
#' gislason(length,params)
#' }
setMethod('gislason', signature(length='FLQuant',params='numeric'),
          function(length,params,a=0.55,b=1.44,c=-1.61,...) { 
            res=gislasonFn(length,params)
            res@units='yr^-1'
            res})
setMethod('gislason', signature(length='FLQuant',params='FLPar'),
          function(length,params,a=0.55,b=1.44,c=-1.61,...){   
            res=gislasonFn(length,params)
            res@units='yr^-1'
            res})

gislasonFn<-function(length,params,a=0.55,b=1.44,c=-1.61) {
  
  # Natural mortality parameters from Model 2, Table 1 gislason 2010
  if (!all(c("m1","m2")%in%dimnames(params)$params)){
    
    m1=FLPar(m1= a*(params["linf"]^b)%*%params["k"], iter=dims(params)$iter)
    m2=FLPar(m2=c                           ,        iter=dims(params)$iter)
    params=rbind(params,m1,m2)
  }
  
  params["m1"]%*%(exp(log(length)%*%params["m2"]))}

#' lorenzen
#'
#' Lorenzen natural mortality relationship estimate M as a function of weight. 
#' M=a*wt^b; 
#' 
#' @param wt  mass at which M is to be predicted
#' @param params an \code{FLPar} with two values; i.e. a equal to M at unit mass 
#' and b a power term; defaults are a=0.3 and b=-0.288
#' @param ... any other arguments
#' 
#' @aliases lorenzen lorenzen-method lorenzen,FLQuant,FLPar-method lorenzen,FLQuant,missing-method lorenzen,FLQuant,numeric-method  lorenzen,numeric,missing-method
#' 
#' @export
#' @docType methods
#' @rdname lorenzen
#' 
#' @seealso \code{\link{gislason}}
#'  
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' mass=FLQuant(c( 1.90, 4.23, 7.47,11.48,16.04,20.96,26.07,31.22,
#'                36.28,41.17,45.83,50.20,54.27,58.03,61.48,64.62),
#'              dimnames=list(age=1:16))
#' lorenzen(mass)
#' }

setMethod('lorenzen', signature(wt='FLQuant',params='FLPar'),
          function(wt,params,...){   
            res=params[1]%*%(wt%^%params[2])
            res@units='yr^-1'
            res})
setMethod('lorenzen', signature(wt='FLQuant',params='missing'),
          function(wt,m1=.3,m2=-0.288,...) { 
            res=lorenzenFn(wt,m1=m1,m2=m2)
            res@units='yr^-1'
            res})
setMethod('lorenzen', signature(wt='FLQuant',params='numeric'),
          function(wt,params,...) { 
            res=params[1]*wt^params[2]
            res@units='yr^-1'
            res})
setMethod('lorenzen', signature(wt='numeric',params='missing'),
          function(wt,m1=.3,m2=-0.288,...) { 
            res=lorenzenFn(wt,m1=m1,m2=m2)
            res})

lorenzenFn<-function(wt,m1=.3,m2=-0.288){
  if ("FLPar"%in%is(m2)) res=wt%^%m2  else res=wt^m2
  if ("FLPar"%in%is(m1)) res=m1%*%res else res=m1*res
  res}

m1<-function(m,wt){
  
  fn<-function(x,wt,ref) sum((lorenzen(wt,m1=x)-m)^2)
  
  optimize(fn, c(0, 100), tol=0.0000001,wt=wt,ref=m)$minimum}



#' @title Natural mortality
#' 
#' @description 
#' Methods to provide estimates of natural mortality based on growth and reproduction parameters
#' 
#' @import FLCore 
#' 
#' @param params \code{FLPar}
#' @param ... any other arguments
#' 
#' @aliases gislasen gislasen-method gislasen,FLQuant,FLPar-method 
#'          lorenzen lorenzen-method lorenzen,FLQuant,FLPar-method
#'          roff roff-method roff,FLPar-method 
#'          rikhter rikhter-method rikhter,FLPar-method
#'          rikhter2 rikhter2-method rikhter2,FLPar-method
#'          griffiths griffiths-method griffiths,FLPar-method
#'          djababli djababli-method djababli,FLPar-method
#'          jensen jensen-method jensen,FLPar-method
#'          jensen2 jensen2-method jensen2,FLPar-method
#'          charnov charnov-method charnov,FLQuant,FLPar-method
#'          petersen petersen-method petersen,FLPar-method petersen,FLQuant,FLPar-method
#'          
#' @return returns an object of \code{FLQuant}
#' 
#' @exportMethod roff rikhter rikhter2 griffiths djababli jensen jensen2 charnov petersen
#' @docType methods
#' @rdname m
#' 
#' @details
#' 
#' Natural Mortality
#' For larger species securing sufficient food to maintain a fast growth rate may entail 
#' exposure to a higher natural mortality @gislason2008does. While many small demersal species 
#' seem to be partly protected against predation by hiding, cryptic behaviour, being flat 
#' or by possessing spines have the lowest rates of natural mortality @griffiths2007natural. 
#' Hence, at a given length individuals belonging to species with a high \deqn{L_{\infty}} may 
#' generally be exposed to a higher M than individuals belonging to species with a low \deqn{L_{\infty}}.
#'
#' \deqn{ log(M) = 0.55-1.61log(L) + 1.44log(L_{\infty}) + log(k)}

#' 
#' Functional forms
#' 
#' Many estimators have been propose for M, based on growth and reproduction,
#' 
#' Age at maturity 
#'  \deqn{M=\frac{1.521}{a_{50}^{0.72}}-0.155}
#'  \deqn{M=\frac{1.65}{a_{50}}}
#' 
#' Growth
#' \deqn{M=1.5k}
#' \deqn{M=1.406W_{\infty}^{-0.096}k^{0.78}}
#' 
#' \deqn{M=1.0661L_{\infty}^{-0.1172}k^{0.5092}}
#' Growth and length at maturity
#' 
#' \deqn{M=3kL_{\infty}\frac{(1-\frac{L_{50}}{L_{\infty}})}{L_{50}}}
#' \deqn{M=\frac{\beta k}{e^{k(a_{50}-t_0)}-1}}
#' 
#' Varing by length, weight or age
#' 
#' @seealso \code{\link{gislason}},  \code{\link{lorenzen}}
#' 
#' @examples
#' \dontrun{
#' params=FLPar(FLPar(linf=120,k=.15,t0=-0.1,l50=60,a=0.0001,b=3))
#' age=FLQuant(1:10,dimnames=list(age=1:10))
#' 
#' 
#' roff(params)
#' rikhter(params)
#' rikhter2(params)
#' griffiths(params)
#' djababli(params)
#' jensen(params)
#' jensen2(params)
#' }
setGeneric('roff',      function(params,...)  standardGeneric('roff'))
setGeneric('rikhter',   function(params,...)  standardGeneric('rikhter'))
setGeneric('rikhter2',  function(params,...)  standardGeneric('rikhter2'))
setGeneric('griffiths', function(params,...)  standardGeneric('griffiths'))
setGeneric('djababli',  function(params,...)  standardGeneric('djababli'))
setGeneric('jensen',    function(params,...)  standardGeneric('jensen'))
setGeneric('jensen2',   function(params,...)  standardGeneric('jensen2'))
setGeneric('charnov',   function(len,params,...)  standardGeneric('charnov'))
setGeneric('petersen',  function(wt,params,...)   standardGeneric('petersen'))
setGeneric('chen',      function(age,params,...)  standardGeneric('chen'))

setMethod('roff', signature(params='FLPar'),function(params,...){
  res=(3*params["k"]%*%params["linf"])*(1.0-params["l50"]%/%params["linf"])%/%params["l50"]
  
  dimnames(res)$params="m"
  res})

setMethod('rikhter', signature(params='FLPar'),function(params,...){
  tm=params["t0"]%+%log(1-params["l50"]/(params["linf"]))%/%(-params["k"])
  res=params["b"]%*%params["k"]%/%(exp(params["k"]%*%(tm%-%params["t0"]))-1)

  dimnames(res)$params="m"
  res})

setMethod('rikhter2', signature(params='FLPar'),function(params,...){
  tm=params["t0"]%+%log(1-params["l50"]/(params["linf"]))%/%(-params["k"])
  res=1.521/tm^0.73-0.155
  
  dimnames(res)$params="m"
  res})

setMethod('griffiths', signature(params='FLPar'),function(params,...){
  winf=params["a"]%*%(params["linf"]^params["b"])
  res=(1.406*(winf^-0.096))%*%params["k"]^0.78
  
  dimnames(res)$params="m"
  res})

setMethod('djababli', signature(params='FLPar'),function(params,...){
  res=(1.066*params["linf"]^-0.1172)%*%params["k"]^0.5092
  
  dimnames(res)$params="m"
  res})

setMethod('jensen', signature(params='FLPar'),function(params,...){
  res=1.5*params["k"]
  
  dimnames(res)$params="m"
  res})

setMethod('jensen2', signature(params='FLPar'),function(params,...){
  tm=params["t0"]%+%log(1-params["l50"]/(params["linf"]))%/%(-params["k"])
  res=1.65/tm
  
  dimnames(res)$params="m"
  res})

setMethod('charnov', signature(len="FLQuant",params='FLPar'),function(len,params,...){
  res=params["k"]%*%(params["linf"]%/%len)^1.5

  res})

setMethod('petersen', signature(wt="FLQuant",params='FLPar'),function(wt,params,...){
  1.28*wt^(-0.25)})

setMethod('chen', signature(age="FLQuant",params='FLPar'),function(age,params,...){ #(age,k,t0=-0.1){
  m =params["k"]/(1-exp(-params["k"]%*%(age%-%params["t0"])))
  
  tm =-(1/params["k"])*log(1-exp(params["k"]*params["t0"]))+params["t0"]
  bit=exp(-params["k"]*(tm-params["t0"]))
  
  a0=1-bit
  a1=params["k"]*bit
  a2=-0.5*params["k"]^2*bit
  age.=age>c(tm)
  m[age.] =params["k"]/(a0+a1*(age[age.]-tm)+a2*(age[age.]-tm)^2)
  
  dimnames(m)$params="m"
  return(m)})   


if (FALSE){
  lh=lhPar(FLPar(linf=100))
  lh=rbind(lh,mPar("gislason"))
  
  mFn("gislason")(wt2len(stock.wt(ple4),lh),lh)
  }
