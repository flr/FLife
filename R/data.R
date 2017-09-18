#' cas
#' 
#' A dataset containing lengths by year
#' 
#' @docType data
#' @rdname cas
#' 
#' @seealso \code{\link{teleost}}
#' 
#' @format An \code{data.frame} object
#' \describe{
#'   \item{year}{ of capture}
#'   \item{len}{ at capture}   
#'   \item{n}{ frequency}   
#' }
#' 
#' @name cas
#' @docType data
#' 
NULL

#' teleost
#'
#' A dataset containing life history parameters for a range of teleost species
#' 
#' @docType data
#' 
#' @seealso \code{\link{cas}}
#' 
#' @format A data frame with 139 rows (i.e. species) and 13 variables:
#' \describe{
#'   \item{species,genus,family,order,class}{} 
#'   \item{linf}{asymptotic length of the Von Bertalanffy growth equation}    
#'   \item{k}{rate at which length reaches linf}       
#'   \item{t0}{adjusts the growth equation for the initial size at birth}      
#'   \item{l50}{length at which 50\% are mature}     
#'   \item{a}{scaling factor of the length weight relationship}       
#'   \item{b}{exponent of the length weight relationship}       
#'   \item{temp}{temperature}    
#'   \item{habit}{e.g. pelagic or demeresal}   
#' }
#' @name teleost
#' 
NULL