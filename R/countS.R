#' @title Generates binomial count data for model Count:S with lambda depending 
#' on a covariate.
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site and depending on the covariate value attached to the site. 
#' Then generates binomial random variables for each site, given the number of animals present, 
#' assuming a constant hazard of detection for each animal. 
#'
#' @param param A vector comprised of parameters \eqn{b_0}{b0}, \eqn{b_1}{b1}, 
#' and the log of the detection hazard \eqn{h}{h} (in that order), where the log of the Poisson rate 
#' lambda is equal to \eqn{b_0+b_1x}{b0+b1*x} and \eqn{x}{x} is the covariate \code{covar}.
#' @param R The number of sites.
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param covar A vector covariate of length \code{R} on which the expected number of 
#' animals in the site depends linearly (assumed to be the same for all occasions).
#' 
#' @return Returns an \code{R} by 1 matrix of count values.
#' 
#' @examples 
#' Rsites=100. #number of sites
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
#' x=rnorm(Rsites,0,1)
#' lambda=exp(b0t+b1t*x)
#' mean(lambda)
#' 
#' # true probabilites
#' 1-exp(-ht*Tsearch)      #P(individual detected)
#' 1-exp(-mean(lambda)) #P(site occupancy)
#' paramt=c(b0t, b1t,ht)
#' 
#' # data for Binary:S
#' cnts=as.matrix(generate.countScov(paramt,R=Rsites,Tmax=Tsearch,covar=x))
#' str(cnts)
#' head(cnts)
#' 
#' @export
generate.countScov=function(param, R, Tmax, covar)
{
  b0=param[1]
  b1=param[2]
  h=param[3]
  lam=exp(b0+b1*covar)
  p=1-exp(-h*Tmax)
  yvec=matrix(0,R,1)
  for(i in 1:R){
    n=rpois(1,lam[i])
    yvec[i]=rbinom(1,n,p)}
  return(yvec)
}

#' @title Evaluates the negeative log-likelihood for model Count:S with lambda depending 
#' on a covariate.
#'
#' @description
#' Evaluates the negeative log-likelihood for model Count:S, assuming that lambda 
#' depends on a covariate \code{covar} that is site-dependent, 
#' given initial parameter estimates and binary data from a multiple-occasion survey.
#'
#' @param param A vector comprised of parameters \eqn{b_0}{b0}, \eqn{b_1}{b1}, 
#' and the log of the detection hazard \eqn{h}{h} (in that order), where the log of the Poisson rate 
#' lambda is equal to \eqn{b_0+b_1x}{b0+b1*x} and \eqn{x}{x} is the coaviate \code{covar}.
#' @param R The number of sites.
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by 1 matrix of binary values, with a 1 representing
#' detection and a 0 representing no detection
#' @param covar A vector covariate of length \code{R} on which the expected number of 
#' animals in the site depends linearly (assumed to be the same for all occasions).
#' 
#' @return Returns the negative log-likelihood function evaluated for model Count:S at 
#' the parameter values passed in \code{param}.
#' 
#' @examples 
#' Rsites=100. #number of sites
#' Tsearch=3 #maximum time
#' 
#' #true parameters
#' b0t =0.660878  #true intercept for log(lambda) 
#' b1t =0.255413 #true slope for log(lambda)
#' ht=0.4620981  #rate parameter
#' 
#' # simulation of the covariate x over R sites 
#' # hence the site-dependent abundance for the true parameters b0t and b1t
#' set.seed(123) # for reprodicibility
#' x=rnorm(Rsites,0,1)
#' lambda=exp(b0t+b1t*x)
#' mean(lambda)
#' 
#' # true probabilites
#' 1-exp(-ht*Tsearch)      #P(individual detected)
#' 1-exp(-mean(lambda)) #P(site occupancy)
#' paramt=c(b0t,b1t,ht)
#' 
#' # data for Binary:S
#' bns=as.matrix(generate.countScov(paramt,R=Rsites,Tmax=Tsearch,covar=x))
#' 
#' init.paramt=c(b0t, b1t, log(ht))
#' nll = nll.countScov(param=init.paramt, R=Rsites, Tmax=Tsearch,dat=bns, covar=x)
#' nll
#' 
#' # optimize
#' fit.countScov=optim(init.paramt,nll.countScov,R=Rsites,Tmax=Tsearch,dat=bns,covar=x)
#' estpar.countScov=fit.countScov$par
#' # compare estimates and true parameters
#' c(estpar.countScov[1:2],exp(estpar.countScov[3]))
#' paramt 
#' 
#' @export
nll.countScov=function(param,R,Tmax,dat,covar)
{
  b0=param[1]
  b1=param[2]
  h=exp(param[3])
  yvec=dat
  lam=exp(b0+b1*covar)
  p=1-exp(-h*Tmax)
  loglik=sum(yvec*log(lam))+sum(yvec*log(p))-
    sum(lam*p)-sum(log(factorial(yvec)))
  return(-loglik)
}
