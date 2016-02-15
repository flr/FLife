#' teleost
#'
#' A dataset containing life history parameters for a range of teleost species
#' 
#'
#' @format A data frame with 139 rows (i.e. species) and 13 variables:
#' \describe{
#'   \item{species}{}
#'   \item{order}{}   
#'   \item{genus}{}   
#'   \item{family}{} 
#'   \item{class}{} 
#'   \item{linf}{asymptotic length of the Von Bertalanffy growth equation}    
#'   \item{k}{rate at which length reaches linf}       
#'   \item{t0}{adjust the growth equation for the initial size}      
#'   \item{l50}{length at which 50% are mature}     
#'   \item{a}{scaling factor of the length weight relationship}       
#'   \item{b}{exponent of the length weight relationship}       
#'   \item{temp}{temperature}    
#'   \item{habit}{e.g. pelagic or demeresal}   
#' }
"teleost"


#' lengths
#'
#' A dataset containing lengths by year
#'
#' @format A data frame with 139 rows (i.e. species) and 13 variables:
#' \describe{
#'   \item{length}{}
#'   \item{year}{}   
#'   \item{data}{}   
#' }
#' @source \url{http://iccat.int/en/accesingdb.htm}
"lengths"

