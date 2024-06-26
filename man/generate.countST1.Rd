% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/countST1.R
\name{generate.countST1}
\alias{generate.countST1}
\title{Generates time to detections for model CountT1:S with lambda constant}
\usage{
generate.countST1(param, R, Tmax)
}
\arguments{
\item{param}{A vector comprised of the Poisson rate lambda, and the detection hazard, h.}

\item{R}{The number of sites.}

\item{Tmax}{The survey duration (assumed to be the same for all sites)}
}
\value{
Returns an \code{R} by 2 matrix with first column being the number of detections 
in the site and second column being the shortest time to detection of any animals in the 
site (with zero representing no detections).
}
\description{
Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
number of animals in each site. 
Then generates exponential random variables for each site, assuming a constant hazard 
of detection for each animal.
}
\examples{
# setting
Rsites=100  #number of sites
Tsearch=3 #maximum time
#true parameters
lamt = 2
ht=0.4620981
paramt=c(lamt,ht)

cnts=as.matrix(generate.countST1(paramt,R=Rsites,Tmax=Tsearch))
str(cnts)
head(cnts)

}
