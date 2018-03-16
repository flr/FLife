---
title: "FLife"
subtitle: "WKLIFE Life History Relationships"
date: "28 septiembre, 2017"
author: "Laurence Kell"
output: rmarkdown:::pdf_document
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{FLife-WKLife}
  %\VignetteEncoding{UTF-8}
bibliography: refs.bib
tags: FLPKG FLR
license: Creative Commons Attribution-ShareAlike 4.0 International
---

# FLife package





```r
library(ggplot2)
library(FLife)
library(plyr)
library(reshape)
```



# Life history parameters

```r
data(wklife)

wklife
```

```
                         name           common          area    stock sex
1             Clupea harengus          Herring   Celtic Seas  her-nis   F
2       Pollachius pollachius          Pollack     North Sea pol-nsea   C
3                 Molva molva             Ling        Widely lin-comb   C
4         Sebastes norvegicus        Rose fish      Northern  smn-con   C
5           Mullus surmuletus       Red mullet   Celtic Seas mut-comb   F
6         Scopthalmus maximus           Turbot     North Sea tur-nsea   F
7            Microstomus kitt       Lemon sole     North Sea lem-nsea   C
8  Lepidorhombus whiffiagonis           Megrim     North Sea meg-4a6a   C
9              Ammodytes spp.         Sandeels     North Sea  san-ns4   C
10      Pleuronectes platessa           Plaice   Celtic Seas ple-celt   F
11       Merlangius merlangus          Whiting   Celtic Seas whg-7e-k   F
12   Melanogrammus aeglefinus          Haddock   Celtic Seas had-iris   C
13        Lophius piscatorius White anglerfish   Celtic Seas ang-78ab   C
14        Lophius piscatorius White anglerfish     North Sea ang-ivvi   C
15                   Nephrops        Shellfish Biscay-Iberia nep-2829   F
         a    b lmax  linf  l50 a50    t0     k
1  0.00480 3.20   NA  33.0 23.0  NA    NA 0.606
2  0.00760 3.07   NA  85.6 47.1  NA    NA 0.190
3  0.00360 3.11   NA 119.0 74.0 7.2    NA 0.140
4  0.01780 2.97   NA  50.2 40.3  NA  0.08 0.110
5  0.00570 3.24   NA  47.5 16.9  NA    NA 0.210
6  0.01490 3.08   NA  66.7 34.2 2.2  0.29 0.320
7  0.01230 2.97   NA  37.0 27.0  NA    NA 0.420
8  0.00220 3.34   NA  54.0 23.0 3.0    NA 0.120
9  0.00490 2.78   NA  24.0 12.0  NA    NA 1.000
10 0.01100 2.96   NA  48.0 22.9  NA    NA 0.230
11 0.01030 2.40   NA  38.0 28.0  NA -1.01 0.380
12 0.01130 2.96   NA  79.9   NA 2.0 -0.36 0.200
13 0.01980 2.90  133 105.6 73.0  NA -0.38 0.180
14 0.02970 2.84   NA 106.0 61.0  NA    NA 0.180
15 0.00056 3.03   NA  65.0 30.0  NA    NA 0.065
```


![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

**Figure 1** Pairwise scatter plots of life history parameters.


# Equilibrium Dynamics

Create an `FLPar`

```r
wkpar=as(wklife[,6:13],"FLPar")
attributes(wkpar)[names(wklife)[1:5]]=wklife[,1:5]
```

Then use life history relationships to estimate missing values  

```r
par=lhPar(wkpar)
```

and then to derive vectors for processses such as natural mortality 
 

```r
library(FLBRP)

eql=lhEql(par)
```


```r
sel<-function(x) 
  catch.sel(x)%/%fapex(catch.sel(x))

ggplot(FLQuants(eql,"m","catch.sel"=sel,"mat","catch.wt"))+
  geom_line(aes(age,data,col=attributes(wkpar)$name[iter]))+
  facet_wrap(~qname,scale="free")+
  scale_x_continuous(limits=c(0,15))+ 
  guides(colour=guide_legend(title="Species",title.position="top"))
```

![plot of chunk vectors](figure/vectors-1.png)

**Figure 2** Vectors of m, selection pattern, maturity and weight-at-age.

and estimate equilibrium dynamics and reference points, e.g. for lemon sole


```r
plot(iter(eql,7))
```

![plot of chunk eql](figure/eql-1.png)

**Figure 3** Equilibrium curves for lemon sole.

# Simulation

Create a forward projection, i.e. an `FLStock` from an equilibrium object


```r
lmsl=as(iter(eql,7),"FLStock")

plot(lmsl)
```

![plot of chunk eql-lmsl](figure/eql-lmsl-1.png)

**Figure 4** Simulate a stock with increasing F

## Software Versions

* R version 3.4.1 (2017-06-30)
* FLCore: 2.6.5
* FLPKG: 
* **Compiled**: Thu Sep 28 17:55:31 2017
* **Git Hash**: 0bee1ca

## Author information

**Laurence KELL**. laurie.kell.es


## Acknowledgements

This vignette and many of the methods documented in it were developed under the MyDas project funded by the Irish exchequer and EMFF 2014-2020. The overall aim of MyDas is to develop and test a range of assessment models and methods to establish Maximum Sustainable Yield (MSY) reference points (or proxy MSY reference points) across the spectrum of data-limited stocks.

# References {#References}


