---
title: "Life History Relationships"
date: "28 septiembre, 2017"
output: pdf_document
vignette: >
  %\VignetteIndexEntry{diags}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
bibliography: refs.bib
github_document:
    mathjax: TRUE
bibliography: ref.bib
tags: FLife FLR
license: Creative Commons Attribution-ShareAlike 4.0 International Public License
---






[](#top)

[Quick Start](#Quick Start)

[Methods](#Methods)

[Simulation](#Simulation)

[Estimation](#Estimation)

[More Information](#More)

[References](#References)


Life history traits include growth rate; age and size at sexual maturity; the temporal pattern or schedule of reproduction; the number, size, and sex ratio of offspring; the distribution of intrinsic or extrinsic mortality rates (e.g., patterns of senescence); and patterns of dormancy and dispersal. These traits contribute directly to age-specific survival and reproductive functions.^[http://www.oxfordbibliographies.com/view/document/obo-9780199830060/obo-9780199830060-0016.xml] The **FLife** package has a variety of methods for modelling life history traits and functional forms for processes for use in fish stock assessment and for conducting Management Strategy Evaluation (MSE). 

These relationships have many uses, for example in age-structured population models, functional relationships for these processes allow the calculation of the population growth rate and have been used to to develop priors in stock assesments and to parameterise ecological models. 

The **FLife** package has methods for modelling functional forms, for simulating equilibrium `FLBRP` and dynamic stock objects `FLStock`.

[Back to Top](#top)

# Quick Start {#QuickStart}

This section provide a quick way to get running and overview of what functions are available, their potential use, and where to seek help. More details are given in later sections.

The simplest way to obtain **FLife** is to install it from the `FLR` repository via the R console:


```r
install.packages("FLife", repos = "http://flr-project.org/R")
```

See help(install.packages) for more details.

After installing the **FLife** package, you need to load it

```r
library(FLife)
```

There is an example teleost dataset used for illustration and as a test dataset, alternatively you can load your own data.


```r
data(teleost)
```

The dataset contains life history parameters for a range of bony fish species and families, i.e. von Bertalanffy growth parameters ($L_{\infty}, k, t_0$), length at 50\% mature ($L_{50}$), and the length weight relationship ($a, b$). 

When loading a new dataset it is always a good idea to run a sanity check e.g.


```r
is(teleost)
```

```
[1] "FLPar"     "array"     "structure" "vector"   
```


The `teleost` object can be used to create `vectors` or other `objects with values by age using **FLife** methods, e.g. to construct a growth curve for hutchen (*Hucho hucho*)

```r
vonB(1:10,teleost[,"Hucho hucho"])
```

```
 [1]  29.0  40.8  51.5  61.1  69.9  77.8  84.9  91.4  97.3 102.6
```

### Plotting

Plotting is done using **ggplot2** which provides a powerful alternative paradigm for creating both simple and complex plots in R using the *Grammar of Graphics* ^[Wilkinson, L. 1999. *The Grammar of Graphics*, Springer. [doi 10.1007/978-3-642-21551-3_13](http://dx.doi.org/10.1007/978-3-642-21551-3_13).] The idea of the grammar is to specify the individual building blocks of a plot and then to combine them to create the desired graphic^[<http://tutorials.iq.harvard.edu/R/Rgraphics/Rgraphics.html>].

The **ggplot** methods expects a `data.frame` for its first argument, `data` (this has been overloaded by **ggplotFL** to also accept FLR objects); then a geometric object `geom` that specifies the actual marks put on to a plot and an aesthetic that is "something you can see" have to be provided. Examples of geometic Objects (geom) include points (geom_point, for scatter plots, dot plots, etc), lines (geom_line, for time series, trend lines, etc) and boxplot (geom_boxplot, for, well, boxplots!). Aesthetic mappings are set with the aes() function and, examples include, position (i.e., on the x and y axes), color ("outside" color), fill ("inside" color), shape (of points), linetype and size. 


```r
age=FLQuant(1:20,dimnames=list(age=1:20))
len=vonB(age,teleost)

ggplot(as.data.frame(len))+
  geom_line(aes(age,data,col=iter))+
  theme(legend.position="none")
```

![Von Bertalanffy growth curves.](figure/unnamed-chunk-1-1.png)

[Back to Top](#top)


# Methods {#Methods}

## Life History Parameters

![Relationship between life history parameters in the teleost dataset.](figure/unnamed-chunk-2-1.png)

### Growth
Consider the von Bertalanffy growth equation

$$ L_t = L_\infty (1 - e^{(-kt-t_0)})$$

where $L_t$ is length at time t, $L_\infty$ the asymptotic maximum length, $k$ the growth coefficient,  and $t_0$ the time at which an individual would, if it possible, be of zero length. 

As $L_\infty$ increases $k$ declines. in other words at a given length a large species will grow faster than a small species. for example @gislason2008coexistence proposed the relationship 

$$k=3.15L_{\infty}^{-0.64}$$

There also appears to be empirical relationship between $t_0$ and $L_{\infty}$ and $k$   i.e.

$$log(-t_0) = -0.3922 - 0.2752 log(L_{\infty}) - 1.038 log(k)$$ 

Therefore for a value of $L_{\infty}$ or even $L_{max}$ the maximum size observered as 
$L_{\infty}=0.95L_{max}$ then all the growth parameters can be recovered.

### Maturity

There is also a relationship between $L_{50}$ the length at which 50% of individuals are mature

$$l_{50}=0.72L_{\infty}^{0.93}$$

and even between the length weight relationship

$$W=aL^b$$

### Natural Mortality

For larger species securing sufficient food to maintain a fast growth rate may entail exposure to a higher natural mortality @gislason2008does. While many small demersal species seem to be partly protected against predation by hiding, cryptic behaviour, being flat or by possessing spines have the lowest rates of natural mortality @griffiths2007natural. Hence, at a given length individuals belonging to species with a high $L_{\infty}$ may generally be exposed to a higher M than individuals belonging to species with a low $L_{\infty}$.

$$ log(M) = 0.55-1.61log(L) + 1.44log(L_{\infty}) + log(k)$$

## Functional forms

In **FLIfe** there are methods for creating growth curves, maturity ogives and natural mortality vectors, selection patterns, and other ogives. All these methods are used to create `FLQuant` objects.

### Growth
gompertz, richards, vonB

![plot of chunk growth](figure/growth-1.png)

### Ogives
dnormal, knife, logistic, sigmoid


```r
dnormal( age,FLPar(a1=4,sl=2,sr=5000))
knife(   age,FLPar(a1=4))
logistic(age,FLPar(a50=4,ato95=1,asym=1.0))
sigmoid( age,FLPar(a50=4,ato95=1))
```


### Natural Mortality

Many estimators have been propose for M, based on growth and reproduction,  see @kenchington2014natural.

### Natural Mortality

Many estimators have been propose for M, based on growth and reproduction,  see @kenchington2014natural.

#### Age at maturity $a_{50}$

Rikhter and Efanov

$$M=\frac{1.521}{a_{50}^{0.72}}-0.155$$
Jensen

$$M=\frac{1.65}{a_{50}}$$

#### Growth

Jensen

$$M=1.5k$$

Griffiths and Harrod

$$M=1.406W_{\infty}^{-0.096}k^{0.78}$$
where $W_{\infty}=\alpha L_{\infty}^{\beta}$

Djabali

$$M=1.0661L_{\infty}^{-0.1172}k^{0.5092}$$

#### Growth and length at maturity $L_{50}$ 

Roff

$$M=3kL_{\infty}\frac{(1-\frac{L_{50}}{L_{\infty}})}{L_{50}}$$

Rikhter and Efanov

$$M=\frac{\beta k}{e^{k(a_{50}-t_0)}-1}$$
where $a_{50}=t_0+\frac{log(1-\frac{L_{50}}{L_{\infty}})}{-k}$

#### Varing by length

Gislason

$$M_L=1.73L^{-1.61}L_{\infty}^{1.44}k$$

Charnov

$$M_L=k\frac{L_{\infty}}{L}^{1.5}$$

#### Varying by weight

Peterson and Wroblewsk

$$M_W=1.28W^{-0.25}$$

Lorenzen

$$M_W=3W^{-0.288}$$

#### Senescence

### Conversions

ages, len2wt, wt2len





Generation of missing life history relationships


```r
par=lhPar(FLPar(linf=100))
par
```

```
An object of class "FLPar"
params
     linf         k        t0         a         b     ato95       a50 
 100.0000    0.1653   -0.1000    0.0003    3.0000    1.0000    4.3600 
     asym        bg        m1        m2        a1        sl        sr 
   1.0000    3.0000  217.3564   -1.6100    4.3600    2.0000 5000.0000 
        s         v 
   0.9000 1000.0000 
units:  cm 
```


There are relationships between the life history parameters and size, growth, maturation, natural mortality and productivity, as seen in the following.


### Simulation
lhPar, lhEql
   
## Function Forms

## Population dynamics
### Ecological
leslie, r

#### life history traits

```
An object of class "FLPar"
iters:  145 

params
               linf                   k                  t0 
45.100000(28.02114)  0.246667( 0.17297) -0.143333( 0.13590) 
                l50                   a                   b 
22.100000(11.71254)  0.011865( 0.00776)  3.010000( 0.15271) 
units:  NA 
```



```

#### Natural Mortality
![plot of chunk m-gislason](figure/m-gislason-1.png)

#### Stock recruitment


#### Fishery 


### Reference points
lopt, loptAge

### Density Dependence
matdd, mdd

### Parameter estination
moment, powh 

### Stationarity
rod

### Random variables
rnoise


### Refetrence points

```r
library(FLBRP)
data(ple4)
refs(ple4)
```

```
An object of class "FLPar"
params
     b.msy   b.virgin     b.f0.1     b.fmax   b.spr.30  b.spr.100 
  1.76e+06   5.25e+06   2.56e+06   1.85e+06   1.89e+06   2.40e+06 
   b.f0.1_    b.fmax_  b.spr.30_ b.spr.100_  b.current      s.msy 
  2.12e+06   1.53e+06   1.56e+06   1.99e+06   3.20e+05   1.58e+06 
  s.virgin     s.f0.1     s.fmax   s.spr.30  s.spr.100    s.f0.1_ 
  5.04e+06   2.34e+06   1.64e+06   1.68e+06   2.19e+06   1.94e+06 
   s.fmax_  s.spr.30_ s.spr.100_  s.current      r.msy   r.virgin 
  1.35e+06   1.39e+06   1.81e+06   2.06e+05   1.05e+06   1.13e+06 
    r.f0.1     r.fmax   r.spr.30  r.spr.100    r.f0.1_    r.fmax_ 
  1.26e+06   1.26e+06   1.26e+06   1.26e+06   1.04e+06   1.04e+06 
 r.spr.30_ r.spr.100_  r.current      f.msy    f.crash     f.f0.1 
  1.04e+06   1.04e+06   8.44e+05   1.15e-01   6.44e-01   8.76e-02 
    f.fmax   f.spr.30    f.f0.1_    f.fmax_  f.spr.30_  f.current 
  1.35e-01   1.32e-01   8.76e-02   1.35e-01   1.32e-01   3.56e-01 
     y.msy     y.f0.1     y.fmax   y.spr.30    y.f0.1_    y.fmax_ 
  1.43e+05   1.63e+05   1.72e+05   1.72e+05   1.35e+05   1.42e+05 
 y.spr.30_ y.spr.100_  y.current          r         rc         rt 
  1.42e+05   1.38e+05   9.60e+04   4.42e-01   9.38e-02   3.86e+00 
units:   
```

# Simulation {#Simulation}

## Simulation of equilibrium values and reference points


```r
library(FLBRP)
eql=lhEql(par)

ggplot(FLQuants(eql,"m","catch.sel","mat","catch.wt"))+
  geom_line(aes(age,data))+
  facet_wrap(~qname,scale="free")+
  scale_x_continuous(limits=c(0,15))
```

![Age-vectors of growthm natural mortality, maturity and selection pattern](figure/unnamed-chunk-9-1.png)

![Equilibrium curves and reference points.](figure/unnamed-chunk-10-1.png)


```
An object of class "FLPar"
params
      r      rc     msy    lopt      sk    spr0  sprmsy 
 0.3943  0.1397 53.6441 63.5204  0.1954  0.1208  0.0263 
units:  NA NA NA NA NA NA NA 
```

Creation of FLBRP objects



## Stock recruitment relationships
![Stock recruitment relationships for a steepness of 0.75 and vigin biomass of 1000](figure/fig3-1.png)

![Production curves, Yield v SSB, for a steepness of 0.75 and vigin biomass of 1000.](figure/fig4-1.png)







## Density Dependence

Modelling density dependence in natural mortality and fecundity.



```r
library(FLBRP)
library(FLife)

data(teleost)
par=teleost[,"Hucho hucho"]
par=lhPar(par)
hutchen=lhEql(par)
 
scale=stock.n(hutchen)[,25]%*%stock.wt(hutchen)
scale=(stock.n(hutchen)%*%stock.wt(hutchen)%-%scale)%/%scale
 
m=mdd(stock.wt(hutchen),par=FLPar(m1=.2,m2=-0.288),scale,k=.5)   

ggplot(as.data.frame(m))+
   geom_line(aes(age,data,col=factor(year)))+
   theme(legend.position="none")+
    scale_x_continuous(limits=c(0,15))
```

![Density Dependence in M](figure/m-density-dependence-1.png)


```r
scale=stock.n(hutchen)[,25]%*%stock.wt(hutchen)
scale=(stock.n(hutchen)%*%stock.wt(hutchen)%-%scale)%/%scale

mat=matdd(ages(scale),par,scale,k=.5)   
 
ggplot(as.data.frame(mat))+
    geom_line(aes(age,data,col=factor(year)))+
    theme(legend.position="none")+
    scale_x_continuous(limits=c(0,15))
```

![Density Dependence in M](figure/Maturity-density-dependence-1.png)






## Noise

Methods to simulate random noise with autocorrelation, e.g. by age or cohort 
![plot of chunk fig9](figure/fig9-1.png)

![plot of chunk fig10](figure/fig10-1.png)









## MSE using empirical HCR
![MSE using empirical HCR](figure/fig14-1.png)

[Back to Top](#top)

# Estimation {#Estimation}

Life history parameters can also be used to estimate quantities of use in stock assessment

@beverton1956review developed a method to estimate life history and population parameters length data. e.g. 

\begin{equation}Z=K\frac{L_{\infty}-\overline{L}}{\overline{L}-L^\prime} \end{equation}
Based on which @powell1979estimation developed a method, extended by @wetherall1987estimating, to estimate growth and mortality parameters. This assumes that the right hand tail of a length frequency distribution was determined by the asymptotic length $L_{\infty}$ and the ratio between Z and the growth rate k. 

The Beverton and Holt methods assumes good estimates for K and $L_{\infty}$, while the Powell-Wetherall method only requires an estimate of K, since $L_{\infty}$ is estimated by the method as well as Z/K.These method therefore provide estimates for each distribution of Z/K, if K is unknown and Z if K is known.  
%As well as assuming that growth follows the von Bertalanffy growth function, it is also assumed that the population is in a steady state with constant exponential mortality, no changes in selection pattern of the fishery and constant recruitment. In the Powell-Wetherall method $L^\prime$ can take any value between the smallest and largest sizes. Equation 1 then provides a series of estimates of Z and since 

\begin{equation}\overline{L}-L^\prime=a+bL^{\prime} \end{equation}
 a and b can be estimated by a regression analysis where 
\begin{equation}b=\frac{-K}{Z+K} \end{equation}
 \begin{equation}a=-bL_{\infty} \end{equation}
Therefore plotting $\overline{L}-L^\prime$ against $L^\prime$ therefore provides an estimate of $L_{\infty}$ and Z/K

Plotting $\overline{L}-L^\prime$ against $L^\prime$ provides an estimate of $L_{\infty}$ and Z/k, since $L_{\infty}=-a/b$ and $Z/k=\frac{-1-b}{b}$. If k is known then it also provides an estimate of Z (\textbf{Figure} \ref{fig:15}).


```
  age   obs    hat    sel
1   1 32356 249252 0.0136
2   2 49911 152624 0.0342
3   3 69038  93457 0.0773
4   4 45627  57226 0.0834
5   5 32732  35041 0.0977
6   6  8910  21457 0.0434
```


![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

### Catch curve analysis

```r
data(ple4)
ctc=as.data.frame(catch.n(ple4))
ctc=ddply(ctc,.(year), with, cc(age=age,n=data))
ctc=ddply(transform(ctc,decade=factor(10*(year%/%10))),.(decade,age),with,data.frame(sel=mean(sel)))
ggplot(ctc)+
  geom_line(aes(age,sel,colour=decade))
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png)

[Back to Top](#top)

# More Information {#More}

* You can submit bug reports, questions or suggestions on `FLife` at the `FLife` issue page ^[<https://github.com/lauriekell/FLife/issues>], or on the *FLR* mailing list.
* Or send a pull request to <https://github.com/lauriekell/FLife/>
* For more information on the FLR Project for Quantitative Fisheries Science in R, visit the FLR webpage ^[<http://flr-project.org>].
* The latest version of `FLife` can always be installed using the `devtools` package, by calling

```r
	library(devtools)
	install_github("lauriekell/FLife")
```
`

## Software Versions

* R version 3.4.1 (2017-06-30)
* FLCore: 2.6.5
* FLPKG: 
* **Compiled**: Thu Sep 28 15:14:14 2017
* **Git Hash**: 3b299a8

## Author information

**Laurence KELL**. laurie.kell.es


## Acknowledgements

This vignette and many of the methods documented in it were developed under the MyDas project funded by the Irish exchequer and EMFF 2014-2020. The overall aim of MyDas is to develop and test a range of assessment models and methods to establish Maximum Sustainable Yield (MSY) reference points (or proxy MSY reference points) across the spectrum of data-limited stocks.

# References {#References}


[Back to Top](#top)
