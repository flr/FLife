
roffFn<-function(params){
  
  (3*params["k"]%*%params["linf"])*(1.0-params["lmat"]%/%params["linf"])%/%params["lmat"]
  }

rikhterEfanovFn<-function(params){
  tm=params["t0"]%+%length*(1-params["lmat"]/(params["linf"]))%/%(-params["k"])
  beta%*%params["k"]%/%(exp(params["k"]%*%(tm%-%params["t0"]))-1)}

rikhterEfanov2Fn<-function(params){
  tm=params["t0"]%+%length*(1-params["lmat"]/(params["linf"]))%/%(-params["k"])
  1.521/tm^0.73-0.155}


griffithsHarrodFn<-function(params){

  winf=params["a"]%*%(length^params["b"])
  (1.406%*%winf^-0.096)%*%params["k"]^-0.78  
  }

djababliFn<-function(params){
  
  (1.066*params["linf"]^-0.1172)%*%params["k"]^0.5092}

jensenFn<-function(params){
  
  1.5*params["k"]}

jensen2Fn<-function(params){
  tm=params["t0"]%+%length*(1-params["lmat"]/(params["linf"]))%/%(-params["k"])
  1.65%/%tm}

charnovFn<-function(params){
  
  params["k"]%*%(params["linf"]%/%length)^1.5}

petersenWroblewskiFn<-function(weight){
  1.28*weight^(-0.25)}

chenWatanabeFn=function(age,params) { #(age,k,t0=-0.1){
  m =params["k"]/(1-exp(-params["k"]*(data-params["t0"])))
  
  tm =-(1/params["k"])*log(1-exp(params["k"]*params["t0"]))+params["t0"]
  bit=exp(-params["k"]*(tm-params["t0"]))
  
  a0=1-bit
  a1=params["k"]*bit
  a2=-0.5*params["k"]^2*bit
  age.=age>c(tm)
  m[age.] =params["k"]/(a0+a1*(age[age.]-tm)+a2*(age[age.]-tm)^2)
  
  return(m)}   


#' Adams M. S., (1999) Population Dynamics and Assessment of Skipjack Tuna (Katsuwonus pelamis) in the
#' Maldives. Ph. D dissertation, Renewable Resources Assessment Group. Centre for Environmental
#' Technology, T. H. Huxley School of Environment, Earth Sciences and Engineering, Imperial College of
#' Science Technology and Medicine, Univ. of London, 304 p.
#' Amarasiri C., Joseph L., (1987) Skipjack tuna (Katsuwonus pelamis) - aspects on the biology and fishery from
#' the western and southern coastal waters of Sri Lanka. IPTP. Coll. Vol. Work. Docs. 2, 1-10.
#' Amorim A. F., Antunes S. A., Arfelli C. A., (1981) Relationships of Katsuwonus pelamis Linnaeus, 1758,
#' caught in the south and southeast of Brazil: Length-weight and gilled/gutted weight. Col.Vol. Sci. Pap.
#' ICCAT, 15 (1): 129-134.
#' Andrade H., Oliveira Campos R. (2002) Allometry coefficient variations of the length–weight relationship of
#' skipjack tuna (Katsuwonus pelamis) caught in the southwest South Atlantic. Fisheries Research 55 (2002)
#' 307–312.
#' Anonymous (2000) Annual report of the Inter-American Tropical Tuna Commission, 1998: 357 pp.
#' Ashida H., Tanabe T., Suzuki N., (2010) Assessment male skipjack tuna spawning activity in the tropical
#' western and central Pacific Ocean. WCPFC-SC6-2010/BI- WP-02, 12p.
#' Bard F.X., Antoine L., (1986) Croissance du listao dans l’Atlantique Est. In: Symons, P.E.K., Miyake, P.M.,
#' Sakagawa, G.T. (Eds.), Proc. ICCAT conference on the international skipjack year program, Madrid, pp.
#' 301–308.
#' Batts B. S., (1972) Sexual maturity, fecundity and sex ratios of: the skipjack tuna, Katsuwonus pelamis, in North
#' Carolina waters. Chesapeake Sci. 13, 193-200.
#' Bayliff W. H., (1988) Growth of skipjack, Katsuwonus pelamis and yellowfin, Thunnus albacares, tunas in the
#' eastern Pacific Ocean, as estimated from tagging data. Inter-Amer. Trop. Tuna Comm. Bull. 19(4), 311-358.
#' Bousquet N., Dortel E., Chassot E., Million J., Eveson J.P., and, Hallier J.-P. (2012). Preliminary assessments of
#' tuna natural mortality rates from a Bayesian Brownie-Petersen model. IOTC-2012-WPTT14-41.
#' Brock V. E., (1954) Some aspects of the biology of the aku (Katsuwonus pelamis) in the Hawaiian Islands. Pac.
#' Sci. 8: 94-104.
#' Brouard F., Grandperrin R., Cillaurren E., (1984) Croissance des jeunes thons jaunes (Thunnus albacares) et des
#' bonites (Katsuwonus pelamis) dans le pacifique tropical occidental. ORSTOM de Port-Vila, Notes et
#' Documents d'Oceanographie 10, 23p.
#' Cayre P. (1981) Maturité sexuelle, fécondité et sex-ratio du listao Katsuwonus pelamis des côtes d’Afrique de
#' l’Ouest (0”-20”N), étudiées à partir des débarquements thoniers (1977-1979) au port de Dakar (Sénégal).
#' ICCAT Col. Vol. Sci. Pap. XV(1), 135-149.
#' Cayré P., Diouf T., Fonteneau A., (1986), Analyse des données de marquages et recaptures de listao
#' (Katsuwonus pelamis) réalisés par le Sénégal et la République du Cap-Vert. In: Symons, P.E.K., Miyake,
#' P.M., Sakagawa, G.T. (Eds.), Proc. ICCAT conference on the international skipjack year program, Madrid,
#' pp. 309–316.
#' Cayré P., Farrugio H., (1986), Biologie de la reproduction du listao (Katsuwonus pelamis) de l’océan Atlantique.
#' In: Symons, P.E.K., Miyake, P.M., Sakagawa, G.T. (Eds.), Proc. ICCAT conference on the international
#' skipjack year program, Madrid, pp. 252–272.
#' Charnov E.L., Gislason H. and Pope J.G. (2013) Evolutionary assembly rules for fish life histories. Fish and
#' Fisheries. 14: 213-224.
#' 194Chi K. S., Yang R.-T., (1973) Age and growth of skipjack tuna in the waters around the southern part of Taiwan.
#' [Summ. in Chin.] Nat. Taiwan Univ. Sci. Rep. Acta Oceanogr. Taiwanica 3, 199-221.
#' Chur V. N., Zharov, V. L., (1983) Determination of age and growth of the skipjack tuna, Katsuwonus pelamis
#' (Scombridae) from the southern part of the Gulf of Guinea. J. Ichthyol. 23(3), 53-67.
#' Chu Tien Vinh. (2000) Study on biology of tuna in the south China Sea, Area IV: Vietnamese Waters. Proc. of
#' the SEAFDEC seminar on fishery resources in the south China Sea, Area IV: Vietnamese waters, p. 146-168.
#' De Bryun P., and Murua, H., (2008) Two simple alternative growth models for skipjack tuna (Katsuwonus
#' pelamis) in the Indian Ocean, as estimated from tagging data. IOTC-2008-WPTDA-09, 12p.
#' Dick E.J., MacCall A. D. (2011) Depletion-Based Stock Reduction Analysis: A catch-based method for
#' determining sustainable yields for data-poor fish stocks. Fisheries Research 110 (2011) 331– 341.
#' Djabali F., Mehailia A., Koudil M. and Brahmi B. (1993) Empirical equations for the estimation of natural
#' mortality in Mediterranean teleosts. Naga, the ICLARM Quarterly 16, 35–37.
#' Fonteneau A., Pallares P. (2005) Tuna natural mortality as a function of their age; The bigeye tuna (Thunnus
#' obesus) case. Col. Vol. Sci. Pap. ICCAT, 57(2): 127-141.
#' Gaertner D., Hallier, J-P., Dortel, E., Chassot, E. and Fonteneau, A. (2011). An update of the Indian Ocean
#' skipjack growth curve parameters with tagging data. Some new evidences on area-specific growth rates.
#' IOTC Thirteenth Working Party on Tropical Tuna. 8p.
#' Gaertner D. (2010) Estimates of Historic Changes in Total Mortality and Selectivity for Eastern Atlantic
#' Skipjack (Katsuwonus pelamis) from Length Composition Data. Aquat. Living Resour. 23: 3-11.
#' Gaertner D., Delgado de Molina A., Ariz J., Pianet R., and Hallier J.-P. (2008) Variability of the growth
#' parameters of the skipjack tuna (Katsuwonus pelamis) among areas in the eastern Atlantic: Analysis from
#' tagging data within a meta-analysis approach. Aquat Living Resour. 21: 349-356.
#' Gaertner D., Hallier J-P., Marsac F., and Motah B. A. (2012) Area-specific growth rates of skipjack in the Indian
#' Ocean using tagging data. Indian Ocean Tuna. Tagging Symposium, Grand Bay, Mauritius, 30th October -
#' 2nd November 2012.
#' Gislason H., Daan N., Rice J.C. and Pope J.G. (2010) Size, growth, temperature and the natural mortality of
#' marine fish. Fish and Fisheries 11, 149–158.
#' Grande M. (2013) The reproductive biology, condition and feeding ecology of the skipjack Katsuwonus pelamis,
#' in the Western Indian Ocean. IOTC-2013-WPTT15-INF09. 229 p. + Annexes.
#' Griffiths D. and Harrod C. (2007) Natural mortality, growth parameters, and environmental temperature in fishes
#' revisited. Can. J. Fish. Aquat. Sci. 64, 249–255.
#' Gulland J.A. (1987) Natural mortality and size. Marine Ecology Progress Series 39, 197–199.
#' Hampton J. (2000) Natural mortality rates in tropical tunas: size really does matter. Can. J. Fish. Aquat. Sci. 57,
#' 1002–1010.
#' Hafiz A. (1987) Skipjack fishery in the Maldives. IPTP. Coll. Vol. Work. Docs. 2, 11-22.
#' Hallier J.P., Gaertner D., (2006), Estimated growth rate of the skipjack tuna (Katsuwonus pelamis) from tagging
#' surveys conducted in the Senegalese area (1996-1999) within a meta-analysis framework Col. Vol. Sc. Pap.
#' ICCAT 59(2): 411-420.
#' Hazin F. H. V., Hazin, H. G., Zagaglia, C. R. Travassos, P., Moacir F. G. Júnior (2001) Analyses des captures de
#' la pêche à la senne réalisées par le B.P. Xixili dans l’océan Atlantique Equatorial. Coll. Vol. Sci. Pap. ICCAT
#' 52(2): 488-498.
#' 195Hewitt D.A., Lambert D.M., Hoenig J.M., Lipcius R.N., Bunnell D.B. and Miller T.J. (2007) Direct and indirect
#' estimates of natural mortality for Chesapeake Bay blue crab. Transactions of the American Fisheries Society
#' 136, 1030–1040.
#' Hoenig J.M. (1983) Empirical use of longevity data to estimate mortality rates. Fishery Bulletin 82, 898–903.
#' Jensen A.L. (1996) Beverton and Holt life history invariants result from optimal trade-off of reproduction and
#' survival. Canadian Journal of Fisheries and Aquatic Sciences 53, 820–822.
#' Joseph J., Calkins T. P., (1969) Population dynamics of the skipjack tuna Katsuwonus pelamis of the Eastern
#' Pacific Ocean. Bull. Inter. Amer. Trop. Tuna Comms. 13(1), 273p.
#' Josse E., Le Guen J. C., Kearney R., Lewis A., Smith A., Marec L., Tomlinson P. K., (1979) Growth of skipjack.
#' Occasional Papers SPC, 11, 83.
#' Kenchington T. K. (2013) Natural mortality estimators for information-limited fisheries. Fish and Fisheries,
#' doi:10.1111/faf.12027.
#' Langley A. and Hampton J. (2008) Stock assessment of skipjack tuna in the western and central Pacific Ocean.
#' WCPFC-SC4-2008/SA-WP-4, 74 p.
#' Lee H.H. and Chang Y.J. (2013) Age-structured Natural Mortality for Pacific Blue Marlin Based on Meta-
#' analysis and an Ad Hoc Model. ISC/13/BILLWG-1/07, 19p.
#' Leroy B. (2013) Preliminary results on skipjack (Katsuwonus pelamis) growth. SCTB13 Working Paper, 13p.
#' Lorenzen K. (1996) The relationship between body weight and natural mortality in juvenile and adult fish: A
#' comparison of natural ecosystems and aquaculture. Journal of Fish Biology 49, 627–647.
#' Marcille J., Stéquert B., (1976) Etude préliminaire de la croissance du listao (Katsuwonus pelamis) dans l'ouest
#' de l'Océan Indien tropical. Cah. ORSTOM (Ser. Oceanogr.) 51(44), 201-206.
#' Maunder M. (2012) Status of skipjack tuna in the eastern Pacific ocean in 2011. IATTC, Scientific advisory
#' committee. SAC-03-07a Skipjack assessment 2011, 29p.
#' Matsumoto W. M., Skillman R. A. and Dizon, A. E. (1984): Synopsis of biological data on skipjack tuna.
#' Katsuwonus pelamis. FAO, Fisheries synopsis no. 136, NOAA Tech. Rep. NMFS Circ. 451, 92 pp.
#' Menezes A. A., Aguiar dos Santos R., Fernandes Lin C., Faulstich Neves L., F., Vianna M. (2010)
#' Caracterização das capturas comerciais do bonito listrado, Katsuwonus pelamis, desembarcado em 2007 no
#' Rio de Janeiro, Brasil. Capa 1, 1. (abstract)
#' Mohan M., Kunhikoya K. K., (1985) Age and growth of Katsuwonus pelamis (Linnaeus) and Thunnus albacares
#' (Bonnaterre) from Minicoy waters, In: E. G. Silas, (ed.) Tuna Fisheries of the Exclusive Economic Zone of
#' India: Biology and Stock Assessment 36, CMFRI 143-148 pp. 12.
#' Pagavino M., Gaertner D., (1995) Ajuste de una curva de crecimiento a frecuencias de tallas de atun listado
#' (Katsuwonus pelamis) pescado en el Mar Caribe suroriental. Col. Vol. Sci. Pap. ICCAT 44(2), 303-309.
#' Pascual M., Iribarne O., (1993). How good are empirical predictions of natural mortality? Fish. Res. 16, 17–24.
#' Pauly D. (1980) On the interrelationships between natural mortality, growth parameters, and mean
#' environmental temperature in 175 fish stocks. Journal du Conseil International pour l’Exploration de la Mer
#' 39,175–192.
#' Peterson I. and Wroblewski J.S. (1984) Mortality rate of fishes in the pelagic ecosystem. Can. J. Fish. Aquat.
#' Sci. 41, 1117–1120.
#' 196Pianet R. (1974) Relations poids-longueurs des listaos (Katsuwonus pelamis) pêchés dans le secteur de Pointe-
#' Noire. Col. Vol. Sci. Pap. ICCAT, 2: 126-133.
#' Rothschild B. J., (1966) preliminary assessment of the yield potential of the skipjack tuna (Katsuwonus pelamis)
#' in the Central Pacific Ocean. Polyp. BU-215-M. 12p.
#' Sibert J. R., Kearney R. E., Lawson T. A., (1983) Variation in growth increments of tagged skipjack Katsuwonus
#' pelamis. Tuna and Billfish Assessment Programme. South Pacific Commission, Noumea, Technical Report
#' No. 10. 43 pp.
#' Sivasubramaniam K., (1985) Tunas and their fishery in the EEZs of Maldives and Sri Lanka.
#' BOBP/WP/31, Bay of Bengal Programme (BOBP), Madras.
#' Stequert B., Ramcharrun B. (1976) La reproduction du listao (Katsuwonus pelamis) dans le bassin ouest de
#' l’océan Indien. Aqrrat. Living. Resour., 9, 235-247.
#' Tanabe T., Kayama S., Ogura M., (2003) An outline of the growth study on skipjack tuna (Katsuwonus pelamis)
#' in the Western Pacific. Doc. IOTC WPTT-03-17, 14 p.
#' Tandog-Edralin D. D., Cortes-Zaragoza E. C., Dalzell P., Pauly D. (1990) Some aspects of the biology of
#' skipjack (Katsuwonus pelamis) in Philippines waters. Asian Marine Biology, 7: 15-29.
#' Thapanand-Chaidee T., Pudprommarat C. (2010) A Bayesian Approach based on MCMC Simulation to Weight-
#' Length Relationship: A Case Study on Skipjack Tuna, Katsuwonus pelamis (Linnaeus, 1758). The 11th
#' National Conference on Statistics and Applied Statistics 2010, Chiang Mai, May 27-28, 2010.9p.
#' Timohina O. I., Romanov E. V., (1996) Characteristics of ovogenesis and some data on maturation and
#' spawning of skipjack tuna, Katsuwonus pelamis (Linnaeus, 1758) from the Western part of the equatorial
#' zone of the Indian Ocean. 12p.
#' Uchiyama J. H., Struhsaker P., (1981) Age and growth of skipjack tuna Katsuwonus pelamis, and yellowfin tuna,
#' Thunnus albacares, as indicated by daily growth increments of sagittae. Fish. Bull. 79, 151-162.
#' Vetter E.F. (1988) Estimation of natural mortality in fish stocks: a review. Fish. Bull., 86, 25–43.
#' Vilela M. J. A., Castello J. P., (1991) Estudio de la edad y del crecimiento del barrilete Katsuwonus pelamis, en
#' la region Sur y Sudeste de Brasil. Frente Maritimo 9, 29-35.
#' Wankowski J. W. J., (1981) Estimated growth of surface schooling skipjack tuna, Katsuwonus pelamis, and
#' yellowfin tuna, Thunnus albacares, from the Papua New Guinea Region. Fish. Bull. 79(3), 517-532.
#' Wild A., Hampton, J., (1994) A review of the biology and fisheries for skipjack tuna, Katsuwonus pelamis, in
#' the Pacific Ocean. FAO Fish. Tech. Pap. 336 (2), 1-51.
#' Yao M., (1981) Growth of skipjack tuna in the western Pacific Ocean. Bull. Tohoku Reg. Fish. Res. Lab. 43, 71-
#' 82.