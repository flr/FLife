
  library(FLife)
  rm(list = ls(all = TRUE))
  ls()
  #character(0)
  
  mort <- c(0.8033223, 0.4219499, 0.2816508, 0.2300376, 0.2110502, 0.2040652, 0.2014955, 0.2005502, 0.2002024, 0.2000745)
  lastyr <- 100
  par <- FLPar(linf=90, a=0.00001, sl=1, sr=2000, a1=4, s=0.65, v=1000, a50 = 2)
  range <- c(min=1, max=10, minfbar = 3, maxfbar=6, plusgroup=10)
  parg <- lhPar(par)
  eql <- lhEql(parg, range=range, m = mort)

  mFun <- function(x, params){ (0.2+ 1.64* exp(-1 * ages(m(x)))) }
  par  <- FLPar(linf=90, a=0.00001, sl=1, sr=2000, a1=4, s=0.65, v=1000, a50 = 2)
  range<- c(min=1, max=10, minfbar = 3, maxfbar=6, plusgroup=10)

  parg <- lhPar(par)
  eql  <- lhEql(parg, range=range, m = mFun)
