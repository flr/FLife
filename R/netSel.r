#Selectivity model from stecf -A.Tidd Sept 2006

# mesh = Meshsize (mm), circ = Circumference (open meshes), panelm = Panel mesh size (mm) 
# distn = (L)Distance from codline (m), twine = Twine diameter (mm) 
# lift = Lifting bag (Yes = 1, No = 0)
# change species (cod,had,whg)

netPar<-function (species="cod",
                  mesh, circ, 
                  panelm=90,distn=12, twine=4,
                  lift=1, a = 9.193 ,b = 0.333,
                  d = -1.206, e = -0.0838, f = -2.677, 
                  g =0.8, i= 0.5, j=0.025)
   {
   #Panel factor 
   #NB: Panel factor only works when:
   #(L)Distance from codline (m)<=16m and Ratio >0.8. Outside these ranges the Panel factor = 1.
   #The twine diameter term only applies when the twine diameter is >= 4mm.
                
   panel          <- g - (distn/60) + i * (mesh/panelm) + j * ((twine*2)-4)   #1 #read NB!!
   ratio          <- mesh/panelm
  
   l50.had        <- a + (b * mesh) + (d * twine) + (e * circ) + (f * lift)
   sr.had         <- (0.0713 * mesh) + (-0.4978 * twine)                   
   l25.had        <- l50.had - sr.had/2         
   l50pan.had     <- panel*l50.had 
   l25pan.had     <- l50pan.had - sr.had/2
   
   #All models based on had l50 equation
   
   res=switch(species,
          "had"={c('l50'=l50.had, 'sr'=sr.had, 'l25'=l25.had,
                   'l50pan'=l50pan.had, 'l25pan'=l25pan.had, 'ratio'=ratio)},
          "cod"={l50.cod        <- l50.had * 1.118
                 sr.cod         <- l50.cod * 0.213
                 l25.cod        <- l50.cod - sr.cod/2
                 l50pan.cod     <- l50.cod*panel
                 l25pan.cod     <- l50pan.cod - sr.cod/2
                 c('l50'=l50.cod, 'sr'=sr.cod, 'l25'=l25.cod,
                   'l50pan'=l50pan.cod, 'l25pan'=l25pan.cod, 'ratio'=ratio)},
          "whg"={l50.whg        <- l50.had * 1.162
                 sr.whg         <- l50.whg * 0.286
                 l25.whg        <- l50.whg - sr.whg/2
                 l50pan.whg     <- l50.whg * panel 
                 l25pan.whg     <- l50pan.whg - sr.whg/2
                 whg            <- c('l50'=l50.whg, 'sr'=sr.whg, 'l25'=l25.whg,
                                     'l50pan.'=l50pan.whg, 'l25pan'=l25pan.whg, 'ratio'=ratio)})
   FLPar(res,units="NA")}

netSel<-function(age,params){
    b=2.197/params["sr"]
    a=-params["l50"]*b
    
    length=vonB(age,params)
    
    exp(a+b*length)/(1+ (exp(a+b*length)))}

netSelAge<-function(len,params){
  b=2.197/params["sr"]
  a=-params["l50"]*b
  
  length=vonB(len=len,params)
  
  exp(a+b*length)/(1+ (exp(a+b*length)))}


if(FALSE){
   data(ple4)
   
   params=FLPar(l50=35,sr=5,t0=-.5,k=0.2,linf=100)
   netSel(ages(mat(ple4)),params)

   neph70  <- netPar("cod", 70,100,90,12)
   neph80  <- netPar("cod", 80,100,90,12)
   neph90  <- netPar("cod", 90,100,90,12)
   whit100 <- netPar("cod",100,100, 0,12)
   whit120 <- netPar("cod",120,100, 0,12)
   flat80  <- netPar("cod", 80,100, 0,12,6)
   flat120 <- netPar("cod",120,100, 0,12,6)
   
   netSel(ages(mat(ple4)),rbind(neph80,FLPar(a=1,b=3,t0=-0.3,k=.3,linf=100)))
   }