#' 
invLogit<-function(y) 0.2+exp(y)%/%(1+exp(y))

steepness<-function(params,a=2.706,b=-3.698){
  
  y=2.706-3.698*params["l50"]%/%params["linf"]

  invLogit(y)}
