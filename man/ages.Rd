% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ages.R
\docType{methods}
\name{ages,FLQuant-method}
\alias{ages,FLQuant-method}
\alias{ages}
\alias{ages-method}
\alias{ages,FLCohort-method}
\title{Fills an \code{FLQuant} with ages}
\usage{
\S4method{ages}{FLQuant}(object)
}
\arguments{
\item{object}{\code{FLQuant} or \code{FLCohort}}

\item{...}{any other arguments}
}
\value{
\code{FLQuant} or \code{FLCohort} depending on what the first argument was
}
\description{
Creates \code{FLQuant} and \code{FLCohort} with ages as entries
}
\examples{
\dontrun{
data(ple4)
ages(m(ple4))}
}
\seealso{
\code{\link{knife}} \code{\link{gascuel}} \code{\link{sigmoid}} \code{\link{gompertz}} \code{\link{vonB}} \code{\link{dnormal}} \code{\link{logistic}}
}
