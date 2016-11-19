library(plyr)
g<-read.csv("/home/laurie/Desktop/flr/generic.txt",sep="'")
mth=sort(unique(as.character(g[,2])))[,-1]

m_ply(mth,function(x){
  
  cat("setGeneric('",x,"',function(object,...) standardGeneric('",x,"'))","\n",
      "#' Sum of vector elements.\n",
      append=TRUE,
      file="/home/laurie/Desktop/flr/FLife/R/generic.RRR")
})
  