% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/countST.R
\name{nll.countST}
\alias{nll.countST}
\title{Evaluates the negative log-likelihood for model CountT:S with lambda constant}
\usage{
nll.countST(param, R, Tmax, dat)
}
\arguments{
\item{param}{A vector comprised of the log of the Poisson rate lambda, and the 
log of the detection hazard, h.}

\item{R}{The number of sites.}

\item{Tmax}{The survey duration (assumed to be the same for all sites)}

\item{dat}{An \code{R} by 2 matrix with first column being the number of detections 
in the site and second column being the **sum** of times to detection for all detected 
animals with in the site (with zero representing no detections).}
}
\value{
Returns the negative log-likelihood function evaluated at the parameter values
passed in \code{param}.
}
\description{
Evaluates the negeative log-likelihood for model CountT:S, assuming constant lambda, 
given initial parameter estimates and count data from a single-occasion survey.
}
\examples{
# setting
Rsites=100  #number of sites
Tsearch=3 #maximum time
# true parameters
lamt = 2
ht=0.4620981
paramt=c(lamt,ht)

# data for Binary:M
set.seed(123) # for reproducibility
cnts=as.matrix(generate.countST(paramt,R=Rsites,Tmax=Tsearch))
# optimize
init.paramt=c(log(lamt),log(ht))
fit.countST=optim(init.paramt,nll.countST,R=Rsites,Tmax=Tsearch,dat=cnts)
estpar.countST=fit.countST$par
# compare estimates and true parameters
exp(estpar.countST)
paramt 

}
