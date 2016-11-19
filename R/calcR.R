## calculates slope at origin of Yield/Stock curve
calcR=function(object,val=0.01){
    gradYS=function(x,rp) {
    
              dimnames(refpts(rp))$refpt[5]="fcrash"
              refpts(rp)@.Data[5,,1]       =c(NA,NA,NA,NA,x,NA,NA,NA)
              grad                         =c(computeRefpts(rp)@.Data["fcrash","yield",1])
              
              return(grad)}
    
     if(is.na(refpts(object)["crash","harvest"])) val=refpts(object)["msy","biomass"]*.01
    
    res=grad(gradYS, val, rp=object)
    
    return(res)}