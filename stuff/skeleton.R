library(FLCore)
library(FLBRP)
library(biodyn)
library(Rcpp)
library(roxygen2)

Rcpp.package.skeleton(name          = "FLife", 
                      list          = character(),
                      environment   = .GlobalEnv, 
                      path          = "/home/laurie/Desktop/flr/git", 
                      force         = TRUE,
                      code_files    = paste("/home/laurie/Desktop/flr/git/lh/R",
                                           c("lh-growth-vonB.R"),sep="/"), 
                      cpp_files     = c("/home/laurie/Desktop/flr/git/biodyn/src/biodyn-logistic.cpp",
                                        "/home/laurie/Desktop/flr/git/biodyn/src/biodyn-logistic.h"),
                      example_code  = TRUE, 
                      attributes    = TRUE, 
                      module        = FALSE,
                      author        = "Laurence Kell",
                      maintainer    = "Laurence Kell",
                      email         = "laurie@seaplusplus.co.uk",
                      license       = "GPL (>= 2)")

