library(plyr)
g<-read.csv("/home/laurie/Desktop/flr/generic.txt",sep="'")
mth=sort(unique(as.character(g[,2])))[,-1]

m_ply(mth,function(x){
  
  cat("setGeneric('",x,"',function(object,...) standardGeneric('",x,"'))","\n",
      "#' Sum of vector elements.\n",
      append=TRUE,
      file="/home/laurie/Desktop/flr/FLife/R/generic.RRR")
})
  
  "#' Sum of vector elements.\n
   #'\n 
   #' \\code{sum} returns the sum of all the values present in its arguments./n
   #'\n 
   #' This is a generic function: methods can be defined for it directly\n
   #' or via the \\code{\\link{Summary}} group generic. For this to work properly,\n
   #' first argument.\n
   #' \n
   #' @param\n
   #' \n
   #' @return\n
   #' the arguments \\code{...} should be unnamed, and dispatch is on the\n
   #' \n
   #' @export\n
   #' \n
   #' @rdname\n
   #' \n
   #' @examples\n
   #' \\dontrun{\n
   #'}",sep="",append=TRUE,
   file="/home/laurie/Desktop/flr/FLife/R/generic.RRR",
   )})