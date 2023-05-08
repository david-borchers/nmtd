#' @title Generates time to detections for model CountT:S with lambda constant
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site. 
#' Then generates exponential random variables for each site, assuming a constant hazard 
#' of detection for each animal. 
#'
#' @param param A vector comprised of the Poisson rate lambda, and the detection hazard, h.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' 
#' @return Returns an \code{R} by \code{J} matrix of counts.
#' 
#' @examples # setting
#' Rsites=100  #number of sites
#' Tsearch=3 #maximum time

#' #true parameters
#' lamt = 2
#' ht=0.4620981
#' paramt=c(lamt,ht)
#' 
# # data for Binary:M
#' cnts=as.matrix(generate.countST(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' str(cnts)
#' head(cnts)
#' 
#' @export
generate.countST=function(param, R, J, Tmax)
{
  lam=param[1]
  h=param[2]
  yt=matrix(0,R,2);
  for(i in 1:R)
  { 
    n=rpois(1,lam)
    if( n > 0)
    {
      tvec=sort(rexp(n,h))
      tvecltT=as.matrix(tvec[tvec < Tmax])
      y=nrow(tvecltT)
      tsum=sum(tvecltT)
      yt[i,]=c(y,tsum)
    }
  }
  return(yt)
}

#' @title Evaluates the negeative log-likelihood for model Count:ST with lambda constant
#'
#' @description
#' Evaluates the negeative log-likelihood for model Count:ST, assuming constant lambda, 
#' given initial parameter estimates and count data from a single-occasion survey.
#'
#' @param param A vector comprised of the log of the Poisson rate lambda, and the 
#' log of the detection hazard, h.
#' @param R The number of sites.
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by 1 matrix of counts.
#' 
#' @return Returns the negative log-likelihood function evaluated at the parameter values
#' passed in \code{param}.
#' 
#' @examples # setting
#' Rsites=100  #number of sites
#' Tsearch=3 #maximum time

#' # true parameters
#' lamt = 2
#' ht=0.4620981
#' paramt=c(lamt,ht)
#' 
#' # data for Binary:M
#' set.seed(123) # for reproducibility
#' cnts=as.matrix(generate.countST(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' # optimize
#' init.paramt=c(log(lamt),log(ht))
#' fit.countST=optim(init.paramt,nll.countST,R=Rsites,J=Jsites,Tmax=Tsearch,dat=cnts)
#' estpar.countST=fit.countST$par
#' # compare estimates and true parameters
#' exp(estpar.countST)
#' paramt 
#' 
#' @export
nll.countST=function(param, R, J, Tmax, dat)
{
  lam=exp(param[1])
  h=exp(param[2])
  loglik=matrix(0,R,1)
  tvec=dat
  for(i in 1:R)
  {y=tvec[i,1]
  tsum=tvec[i,2]
  if(y == 0)
  {loglik[i]=loglik[i]-lam*(1-exp(-h*Tmax))}
  else
  {loglik[i]=loglik[i]+y*log(h*lam)-lam*(1-exp(-h*Tmax))-h*tsum}
  }
  return(-sum(loglik))
}

#' @title Generates time to detection for model Count:ST with lambda depending 
#' on a covariate.
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site and depending on the covariate value attached to the site. 
#' Then generates exponential random variables for each anima for each of \code{J} 
#' occasions that each site is surveyed, assuming a constant hazard of detection for 
#' each animal. 
#'
#' @param param A vector comprised of parameters \eqn{b_0}{b0}, \eqn{b_1}{b1}, 
#' and the log of the detection hazard (in that order), where the log of the Poisson rate 
#' lambda is equal to \eqn{b_0+b_1x}{b0+b1*x} and \eqn{x}{x} is the coaviate \code{covar}.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param covar A vector covariate of length \code{R} on which the expected number of 
#' animals in the site depends linearly (assumed to be the same for all occasions).
#' 
#' @return Returns an \code{R} by \code{J} matrix of counts.
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
#' # simulation of the covariate x over R sites 
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
#' # data for Binary:M
#' cnts=as.matrix(generate.countSTcov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' str(cnts)
#' head(cnts)
#' 
#' @export
generate.countSTcov=function(param,R,J,Tmax,covar)
{
  b0=param[1]
  b1=param[2]
  h=param[3]
  lam=exp(b0+b1*covar)
  yt=matrix(0,R,2);
  for(i in 1:R)
  { n=rpois(1,lam[i])
  if( n > 0)
  {
    tvec=sort(rexp(n,h))
    tvecltT=as.matrix(tvec[tvec < Tmax])
    y=nrow(tvecltT)
    tsum=sum(tvecltT)
    yt[i,]=c(y,tsum)
  }
  }
  return(yt)
}

#' @title Evaluates the negeative log-likelihood for model Count:M with lambda depending 
#' on a covariate.
#'
#' @description
#' Evaluates the negeative log-likelihood for model Count:M, assuming that lambda 
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
#' @param covar A vector covariate of length \code{R} on which the expected number of 
#' animals in the site depends linearly (assumed to be the same for all occasions).
#' 
#' @return Returns the negative log-likelihood function for model Count:M evaluated at 
#' the parameter values passed in \code{param}.
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
#' # simulation of the covariate x over R sites 
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
#' # data for Count:ST
#' cntm=as.matrix(generate.countSTcov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' 
#' init.paramt=c(b0t, b1t, log(ht))
#' nll = nll.countSTcov(param=init.paramt, R=Rsites, J=Jsites, Tmax=Tsearch,dat=cntm, covar=x)
#' nll
#' 
#' # optimize
#' fit.countSTcov=optim(init.paramt,nll.countSTcov,R=Rsites,J=Jsites,Tmax=Tsearch,dat=cntm,covar=x)
#' estpar.countSTcov=fit.countSTcov$par
#' # compare estimates and true parameters
#' c(estpar.countSTcov[1:2],exp(estpar.countSTcov[3]))
#' paramt 
#' 
#' @export
nll.countSTcov=function(param,R,J,Tmax,dat,covar)
{
  b0=param[1]
  b1=param[2]
  h=exp(param[3])
  lam=exp(b0+b1*covar)
  loglik=matrix(0,R,1)
  tvec=dat
  for(i in 1:R)
  {
    y=tvec[i,1]
    tsum=tvec[i,2]
    if(y == 0)
    {loglik[i]=loglik[i]-lam[i]*(1-exp(-h*Tmax))}
    else
    {loglik[i]=loglik[i]+y*log(h*lam[i])-lam[i]*(1-exp(-h*Tmax))-h*tsum}
  }
  return(-sum(loglik))
}
