% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/countM.R
\name{generate.countM}
\alias{generate.countM}
\title{Generates count data for model CountT:S with lambda constant}
\usage{
generate.countM(param, R, J, Tmax)
}
\arguments{
\item{param}{A vector comprised of the Poisson rate lambda, and the detection hazard, h.}

\item{R}{The number of sites.}

\item{J}{The number of occasions (assumed the same for all sites).}

\item{Tmax}{The survey duration (assumed to be the same for all sites)}
}
\value{
Returns an \code{R} by \code{J} matrix of counts.
}
\description{
Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
number of animals in each site. 
Then generates binomial random variables for each 
site on each of \code{J} occasions, given the number of animals present, 
assuming a constant hazard of detection for each animal.
}
\examples{
# setting
Rsites=100  #number of sites
Jsites=5    #multiple visits
Tsearch=3 #maximum time
#true parameters
lamt = 2
ht=0.4620981
paramt=c(lamt,ht)

# data for Count:M
cntm=as.matrix(generate.countM(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
str(cntm)
head(cntm)

}
