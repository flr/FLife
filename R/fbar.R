setGeneric("fbar<-", function(object, value){
  standardGeneric("fbar<-")})

# fbar {{{
setMethod("fbar<-", signature(object="FLBRP", value="numeric"),
	function(object, value) {

		fbar(object) <- FLQuant(value, quant="age")

		return(object)
  }
)

setMethod("fbar<-", signature(object="FLBRP", value="FLQuant"),
	function(object, value) {

		fbar(object)<-value

		return(object)
  }
) # }}}
