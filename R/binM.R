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


#' @title Generates binary data for model Binary:M with lambda depending 
#' on a covariate.
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site and depending on the covariate value attached to the site. 
#' Then generates Bernoulli random variables for each of \code{J} occasions that each site 
#' is surveyed, assuming a constant hazard of detection for each animal. 
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
#' paramt=c(b0t, b1t,ht)
#' 
#' # data for Binary:M
#' bnm=as.matrix(generate.binMcov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' str(bnm)
#' head(bnm)
#' 
#' @export
generate.binMcov=function(param,R,J,Tmax,covar)
{
  b0=param[1]
  b1=param[2]
  h=param[3]
  lambda=exp(b0+b1*covar)
  ymat=matrix(0,R,J)
  for(i in 1:R)
  {
    n = rpois(1,lambda[i])
    p=1-exp(-h*n*Tmax)
    for(j in 1:J) {ymat[i,j]=rbinom(1,1,p)}
  }
  return(ymat)
}


#' @title Evaluates the negeative log-likelihood for model Binary:M with lambda depending 
#' on a covariate.
#'
#' @description
#' Evaluates the negeative log-likelihood for model Binary:M, assuming that lambda 
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
#' @return Returns the negative log-likelihood function for model Binary:M evaluated at 
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
#' x=rnorm(Rsites,0,1)
#' lambda=exp(b0t+b1t*x)
#' mean(lambda)
#' 
#' # true probabilites
#' 1-exp(-ht*Tsearch)      #P(individual detected)
#' 1-exp(-mean(lambda)) #P(site occupancy)
#' paramt=c(b0t,b1t,ht)
#' 
#' # data for Binary:M
#' bnm=as.matrix(generate.binMcov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' 
#' init.paramt=c(b0t, b1t, log(ht))
#' nll = nll.binMcov(param=init.paramt, R=Rsites, J=Jsites, Tmax=Tsearch,dat=bns, covar=x)
#' nll
#' 
#' # optimize
#' fit.binMcov=optim(init.paramt,nll.binMcov,R=Rsites,J=Jsites,Tmax=Tsearch,dat=bns,covar=x)
#' estpar.binMcov=fit.binMcov$par
#' # compare estimates and true parameters
#' c(estpar.binMcov[1:2],exp(estpar.binMcov[3]))
#' paramt 
#' 
#' @export
nll.binMcov=function(param,R,J,Tmax,dat,covar)
{
  b0=param[1]
  b1=param[2]
  h=exp(param[3])
  J1=rowSums(dat)
  J0=J-J1
  lam=exp(b0+b1*covar)
  loglik=matrix(0,R,1)
  # note that if J1[i]=0 the expresson reduces correctly to
  # exp(-lam[1-e^(-hTmax J)])
  for(i in 1:R)
  { 
    J1i=J1[i];
    J0i=J0[i];
    term=0;
    for(j in 0:J1i)
    {
      a=exp(-h*Tmax*(j+J0i));
      term=term+(-1)^j*choose(J1i,j)*(exp(lam[i]*a))
    }
    loglik[i]=-lam[i]+log(term)+log(choose(J,J1i))
  }
return(-sum(loglik))
}