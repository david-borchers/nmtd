% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PcountS.R
\name{generate.PcountS}
\alias{generate.PcountS}
\title{Generates count data for model PCount:S with lambda constant}
\usage{
generate.PcountS(param, R, Tmax)
}
\arguments{
\item{param}{A vector comprised of the Poisson rate lambda, and the detection hazard, h.}

\item{R}{The number of sites.}

\item{Tmax}{The survey duration (assumed to be the same for all sites)}
}
\value{
Returns an \code{R} by 1 matrix of counts.
}
\description{
Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
number of animals in each site. Then generates Poisson random variables for each 
site, given the number of animals present, assuming a constant hazard of detection for 
each animal.
}
\examples{
# seed
set.seed(123)
# setting
R=100. #number of sites
Tsearch=3 #maximum time
#true parameters
lamt=2
gt=0.4620981
paramt=c(lamt,gt)

ymat=generate.PcountS(paramt, R, Tmax)
str(ymat)
head(ymat)

}
