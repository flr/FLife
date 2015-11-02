loptFn=function(x,params,age=0:200,
                par){
    
  len   =vonB(params,FLQuant(age,dimnames=list(age=age)))
  m     =fnM(params,len)
  mCum  =FLQuant(aaply(m,2,cumsum),dimnames=dimnames(m))
  
  a =max(0,invVonBFn(params,x))

  a_=floor(a)
  
  n =exp(-mCum[ac(a_)]-m[ac(a_)]*(a-a_))
  c(n*len2wt(params,x))}

#params=par[,10]
#params["t0"]=-.1
#loptFn(x,params)

lopt=function(params,
              fnM=function(par,len) 
                    0.55*(len^-1.66)%*%(par["linf"]^1.44)%*%par["k"])
    optimise(loptFn,c(.01,c(params["linf"])*.99),par=params,maximum=TRUE,fnM=fnM)$maximum

loptFn2<-function(params) params["linf"]*3/(3+exp(params["m2"])/params["k"])

lopt2=function(params,M) { 
  if ("FLPar"   %in% is(M)) return((params["linf"]*3)%/%(3+(M%/%params["k"])))
  if ("numeric" %in% is(M)) return((params["linf"]*3)%/%(3+(M/params["k"])))}

lopt3=function(x){
  
  fbar(x)          =fbar(x)[,1]*0
  stock.n(x)%*%stock.wt(x)}

  