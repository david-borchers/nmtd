#' @title Generates binary data for model Binary:S with lambda depending 
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
#' @param covar A vector covariate of length \code{R} on which the hazard depends 
#' linearly (assumed to be the same for all occasions).
#' 
#' @return Returns an \code{R} by \code{J} matrix of binary values, with a 1 representing
#' detection and a 0 representing no detection.
#' 
#' @examples 
#' Rsites=100. #number of sites
#' Jsites=5    #number of visits
#' Tsearch=3 #maximum time
#' 
#' #true parameters
#' b0t =0.660878  #true intercept for log(lambda) 
#' b1t =0.255413 #true slope for log(lambda)
#' ht=0.4620981  #rate parameter
#' 
#' # simulation of the covariate vegHt over R sites 
#' # hence the site-dependent abundance for the true parameters b0t and b1t
#' set.seed(123) # for reprodicibility
#' x=rnorm(R,0,1)
#' lambda=exp(b0t+b1t*x)
#' mean(lambda)
#' 
#' # true probabilites
#' 1-exp(-ht*Tmax)      #P(individual detected)
#' 1-exp(-mean(lambda)) #P(site occupancy)
#' paramt=c(b0t, b1t,ht)
#' 
#' # data for Binary:S
#' bns=as.matrix(generate.binS(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' str(bns)
#' head(bns)
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
#' Evaluates the negeative log-likelihood for model Binary:S, assuming that lambda 
#' depends on a covariate \code{covar} that is site-dependent, 
#' given initial parameter estimates and binary data from a multiple-occasion survey.
#'
#' @param param A vector comprised of parameters \eqn{b_0}{b0}, \eqn{b_1}{b1}, 
#' and the log of the detection hazard (in that order), where the log of the Poisson rate 
#' lambda is equal to \eqn{b_0+b_1x}{b0+b1*x} and \eqn{x}{x} is the coaviate \code{covar}.
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
#' Rsites=100. #number of sites
#' Jsites=5    #number of visits
#' Tsearch=3 #maximum time
#' 
#' #true parameters
#' b0t =0.660878  #true intercept for log(lambda) 
#' b1t =0.255413 #true slope for log(lambda)
#' ht=0.4620981  #rate parameter
#' 
#' # simulation of the covariate vegHt over R sites 
#' # hence the site-dependent abundance for the true parameters b0t and b1t
#' set.seed(123) # for reprodicibility
#' x=rnorm(R,0,1)
#' lambda=exp(b0t+b1t*x)
#' mean(lambda)
#' 
#' # true probabilites
#' 1-exp(-ht*Tmax)      #P(individual detected)
#' 1-exp(-mean(lambda)) #P(site occupancy)
#' paramt=c(b0t,b1t,ht)
#' 
#' # data for Binary:S
#' bns=as.matrix(generate.binS(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' 
#' init.paramt=c(b0t, b1t, log(ht))
#' nll = nll.binS(param=init.paramt, R=Rsites, J=Jsites, Tmax=Tsearch,dat=bns, covar=x)
#' nll
#' 
#' # optimize
#' fit.binS=optim(init.paramt,nll.binS,R=Rsites,J=Jsites,Tmax=Tsearch,dat=bns,covar=x)
#' estpar.binS=fit.binS$par
#' # compare estimates and true parameters
#' c(estpar.binS[1:2],exp(estpar.binS[3]))
#' paramt 
#' 
#' @export
nll.binS=function(param,R,J,Tmax,dat,covar)
{
  b0=param[1]
  b1=param[2]
  h=exp(param[3])
  lam=as.matrix(exp(b0+b1*covar))
  yvec=dat
  loglik=sum(-(1-yvec)*lam*(1-exp(-h*Tmax))+
               yvec*log(1-exp(-lam*(1-exp(-h*Tmax)))) )
  return(-loglik)
}