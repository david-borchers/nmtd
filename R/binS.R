#' @title Generates binary data for model Binary:M with lambda depending 
#' on a scalar covariate.
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
#' @param covar A (scalar) covariate on which the hazard depends linearly (assumed to be 
#' the same for all sites).
#' 
#' @return Returns an \code{R} by \code{J} matrix of binary values, with a 1 representing
#' detection and a 0 representing no detection.
#' 
#' @examples 
#' ###### NOT YET UPDATED #########
#' # setting
#' Rsites=100  #number of sites
#' Jsites=5    #multiple visits
#' Tsearch=3 #maximum time

#' #true parameters
#' lamt = 2
#' ht=0.4620981
#' paramt=c(lamt,ht)
#' 
# data for Binary:M
#' bnm=as.matrix(generate.binM(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' str(bnm)
#' head(bnm)
#' 
#' @export
generate.binS=function(param,R,J,Tmax,covar)
{
  b0=param[1]
  b1=param[2]
  h=param[3]
  lambda=exp(b0+b1*covar)
  yvec=matrix(0,R,1)
  n = rpois(R, lambda)
  pvec=1-exp(-h*n*Tmax)
  for(i in 1:R)
  {yvec[i]=rbinom(1,1,pvec[i])}
  data=cbind(yvec)
  return(data)
}

#' @title Evaluates the negeative log-likelihood for model Binary:S with lambda depending 
#' on a scalar covariate.
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site. Then generates Bernoulli random variables for each of 
#' \code{J} occasions that each site is surveyed, assuming a constant hazard of detection 
#' for each animal. 
#'
#' @param param A vector comprised of parameters \eqn{b_0}{b0}, \eqn{b_1}{b1}, and the log
#' of the detection hazard (in that order), where the log of the Poisson rate lambda is equal to 
#' \eqn{b_0+b_1x}{b0+b1*x} and \eqn{x}{x} is the coaviate \code{covar}.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by \code{J} matrix of binary values, with a 1 representing
#' detection and a 0 representing no detection
#' @param covar A (scalar) covariate on which the hazard depends linearly (assumed to be 
#' the same for all sites).
#' 
#' @return Returns the negative log-likelihood function evaluated at the parameter values
#' passed in \code{param}.
#' 
#' @examples 
#' ###### NOT YET UPDATED #########
#' # setting
#' Rsites=100  #number of sites
#' Jsites=5    #multiple visits
#' Tsearch=3 #maximum time

#' #true parameters
#' lamt = 2
#' ht=0.4620981
#' paramt=c(lamt,ht)
#' param0=c(log(lamt),log(ht))
#' 
# data for Binary:M
#' bnm=as.matrix(generate.binM(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' # optimize
#' finopt=optim(param0,nll.binM,R=Rsites,J=Jsites,Tmax=Tsearch,dat=bnm)
#' finpar=finopt$par
#' # save
#' exp(finpar[1])
#' exp(finpar[2])
#' 
#' @export
nll.binS=function(param,R,J,Tmax,dat,covar)
{
  b0=param[1]
  b1=param[2]
  h=exp(param[3])
  lam=as.matrix(exp(b0+b1*covar))
  yvec=gdatx
  loglik=sum(-(1-yvec)*lam*(1-exp(-h*Tmax))+
               yvec*log(1-exp(-lam*(1-exp(-h*Tmax)))) )
  return(-loglik)
}