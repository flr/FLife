
setGeneric("production", function(object, ...)
  standardGeneric("production"))

setMethod("production", signature(object="FLBRP"),
          function(object,what="ssb",...){
            
            switch(tolower(substr(what,1,1)),
                   s=catch.obs(object)+window(ssb.obs(object),start=dims(ssb.obs(object))$minyear+1,end=dims(ssb.obs(object))$maxyear+1)-ssb.obs(object),
                   b=catch.obs(object)+window(stock.obs(object),start=dims(stock.obs(object))$minyear+1,end=dims(stock.obs(object))$maxyear+1)-stock.obs(object))
          })

setMethod("production", signature(object="FLStock"),
          function(object,what="ssb",...){
            
            switch(tolower(substr(what,1,1)),
                   s=computeCatch(object)+window(ssb(object),start=dims(ssb(object))$minyear+1,end=dims(ssb(object))$maxyear+1)-ssb(object),
                   b=computeCatch(object)+window(stock(object),start=dims(stock(object))$minyear+1,end=dims(stock(object))$maxyear+1)-stock(object),
                   e={biomass=apply(catch.wt(object)*stock.n(object)*catch.sel(object)%/%apply(catch.sel(object),1,max),seq(6)[-1],sum)
                   computeCatch(object)+window(biomass(object),start=dims(biomass(object))$minyear+1,end=dims(biomass(object))$maxyear+1)-biomass(object)})
          })