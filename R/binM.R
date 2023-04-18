#' @title Generates binary data for model Binary:M with lambda constant
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site. Then generates Bernoulli random variables for each of 
#' \code{J} occasions that each site is surveyed, assuming a constant hazard of detection 
#' for each animal. 
#'
#' @param param A vector comprised of the Poisson rate lambda, and the detection hazard, h.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' 
#' @return Returns an \code{R} by \code{J} matrix of binary values, with a 1 representing
#' detection and a 0 representing no detection.
#' 
#' @examples # setting
#' Rsites=100  #number of sites
#' Jsites=5    #multiple visits
#' Tsearch=3 #maximum time

#' #true parameters
#' lamt = 2
#' ht=0.4620981
#' paramt=c(lamt,ht)
#' 
# # data for Binary:M
#' bnm=as.matrix(generate.binM(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' str(bnm)
#' head(bnm)
#' 
#' @export
generate.binM=function(param,R,J,Tmax)
{
  lam=param[1]
  h=param[2]
  ymat=matrix(0,R,J)
  for(i in 1:R)
  {n = rpois(1, lam)
  p=1-exp(-h*n*Tmax)
  ymat[i,]=rbinom(J,1,p)
  }
  return(ymat)
}



#' @title Evaluates the negeative log-likelihood for model Binary:M with lambda constant
#'
#' @description
#' Evaluates the negeative log-likelihood for model Binary:M, assuming constant lambda, 
#' given initial parameter estimates and binary data from a multiple-occasion survey.
#'
#' @param param A vector comprised of the log of the Poisson rate lambda, and the 
#' log of the detection hazard, h.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by \code{J} matrix of binary values, with a 1 representing
#' detection and a 0 representing no detection
#' 
#' @return Returns the negative log-likelihood function evaluated at the parameter values
#' passed in \code{param}.
#' 
#' @examples # setting
#' Rsites=100  #number of sites
#' Jsites=5    #multiple visits
#' Tsearch=3 #maximum time

#' # true parameters
#' lamt = 2
#' ht=0.4620981
#' paramt=c(lamt,ht)
#' 
#' # data for Binary:M
#' set.seed(123) # for reproducibility
#' bnm=as.matrix(generate.binM(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' # optimize
#' init.paramt=c(log(lamt),log(ht))
#' fit.binM=optim(init.paramt,nll.binM,R=Rsites,J=Jsites,Tmax=Tsearch,dat=bnm)
#' estpar.binM=fit.binM$par
#' # compare estimates and true parameters
#' exp(estpar.binM)
#' paramt 
#' 
#' @export
nll.binM=function(param,R,J,Tmax,dat)
{
  lam=exp(param[1])
  h=exp(param[2])
  J1=rowSums(dat) # number of detections at the R sites
  J0=J-J1
  loglik=matrix(0,R,1)
  for(i in 1:R)
  {
    J1i=J1[i];
    J0i=J0[i];
    term=0;
    for(j in 0:J1i)
    {a=exp(-h*Tmax*(j+J0i));
    term=term+(-1)^j*choose(J1i,j)*(exp(lam*a))
    }
    loglik[i]=-lam+log(term)+log(choose(J,J1i))
  }
  return(-sum(loglik))
}