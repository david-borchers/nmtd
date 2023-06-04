#' @title Generates time to detections for model CountT1:S with lambda constant
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site. 
#' Then generates exponential random variables for each site, assuming a constant hazard 
#' of detection for each animal. 
#'
#' @param param A vector comprised of the Poisson rate lambda, and the detection hazard, h.
#' @param R The number of sites.
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' 
#' @return Returns an \code{R} by 2 matrix with first column being the number of detections 
#' in the site and second column being the shortest time to detection of any animals in the 
#' site (with zero representing no detections).
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
#' cnts=as.matrix(generate.countST1(paramt,R=Rsites,Tmax=Tsearch))
#' str(cnts)
#' head(cnts)
#' 
#' @export
generate.countST1=function(param, R, Tmax)
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
      if(y >0)
      {t1=min(tvecltT)
      yt[i,]=c(y,t1)}
    }
  }
  return(yt)
}



#' @title Evaluates the negative log-likelihood for model CountT1:S with lambda constant
#'
#' @description
#' Evaluates the negeative log-likelihood for model CountT1:S, assuming constant lambda, 
#' given initial parameter estimates and count data from a single-occasion survey.
#'
#' @param param A vector comprised of the log of the Poisson rate lambda, and the 
#' log of the detection hazard, h.
#' @param R The number of sites.
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by 2 matrix with first column being the number of detections 
#' in the site and second column being the time to first detection of an animal 
#' in the site (with zero representing no detections).
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
#' cnts=as.matrix(generate.countST1(paramt,R=Rsites,Tmax=Tsearch))
#' # optimize
#' init.paramt=c(log(lamt),log(ht))
#' fit.countST1=optim(init.paramt,nll.countST1,R=Rsites,Tmax=Tsearch,dat=cnts)
#' estpar.countST1=fit.countST1$par
#' # compare estimates and true parameters
#' exp(estpar.countST1)
#' paramt 
#' 
#' @export
nll.countST1=function(param, R, Tmax, dat)
{
  lam=exp(param[1])
  h=exp(param[2])
  ytmat=dat
  loglik=0
  for(i in 1:R)
  {
    y=ytmat[i,1]
    t1=ytmat[i,2]
    if(y == 0)
    {loglik=loglik-lam*(1-exp(-h*Tmax))}
    else
    {term1=log(h)+y*log(lam)
    term2=(y-1)*log((exp(-h*t1)-exp(-h*Tmax)))
    term3=-lam*(1-exp(-h*Tmax))-h*t1
    term4=-log(factorial(y-1))
    loglik=loglik+term1+term2+term3+term4
    }
  }
  return(-loglik)
}



#' @title Generates time to detection for model CountT1:S with lambda depending 
#' on a covariate.
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site and depending on the covariate value attached to the site. 
#' Then generates exponential random variables for each anima1, assuming a constant hazard of 
#' detection for each animal. 
#'
#' @param param A vector comprised of parameters \eqn{b_0}{b0}, \eqn{b_1}{b1}, 
#' and the log of the detection hazard (in that order), where the log of the Poisson rate 
#' lambda is equal to \eqn{b_0+b_1x}{b0+b1*x} and \eqn{x}{x} is the coaviate \code{covar}.
#' @param R The number of sites.
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param covar A vector covariate of length \code{R} on which the expected number of 
#' animals in the site depends linearly (assumed to be the same for all occasions).
#' 
#' @return Returns an \code{R} by 2 matrix with first column being the number of detections 
#' in the site and second column being the shortest time to detection of any animals 
#' in the site (with zero representing no detections).
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
#' paramt=c(b0t, b1t,ht)
#' 
#' # data for Binary:M
#' cnts=as.matrix(generate.countST1cov(paramt,R=Rsites,Tmax=Tsearch,covar=x))
#' str(cnts)
#' head(cnts)
#' 
#' @export
generate.countST1cov=function(param,R,Tmax,covar)
{
  b0=param[1]
  b1=param[2]
  h=param[3]
  lam=exp(b0+b1*covar)
  yt=matrix(0,R,2);
  for(i in 1:R)
  {
    n=rpois(1,lam[i])
    if( n > 0)
    {
      tvec=sort(rexp(n,h))
      tvecltT=as.matrix(tvec[tvec < Tmax])
      y=nrow(tvecltT)
      if(y >0)
      {
        t1=min(tvecltT)
        yt[i,]=c(y,t1)
      }
    }
  }
  return(yt)
}


#' @title Evaluates the negeative log-likelihood for model CountT1:S with lambda depending 
#' on a covariate.
#'
#' @description
#' Evaluates the negeative log-likelihood for model CountT1:S, assuming that lambda 
#' depends on a covariate \code{covar} that is site-dependent, 
#' given initial parameter estimates and binary data from a multiple-occasion survey.
#'
#' @param param A vector comprised of parameters \eqn{b_0}{b0}, \eqn{b_1}{b1}, 
#' and the log of the detection hazard (in that order), where the log of the Poisson rate 
#' lambda is equal to \eqn{b_0+b_1x}{b0+b1*x} and \eqn{x}{x} is the coaviate \code{covar}.
#' @param R The number of sites.
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by 2 matrix with first column being the number of detections 
#' in the site and second column being the time to first detection of an animal 
#' in the site (with zero representing no detections).
#' @param covar A vector covariate of length \code{R} on which the expected number of 
#' animals in the site depends linearly (assumed to be the same for all occasions).
#' 
#' @return Returns the negative log-likelihood function for model CountT1:S evaluated at 
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
#' # data for CountT1:S
#' cntm=as.matrix(generate.countST1cov(paramt,R=Rsites,Tmax=Tsearch,covar=x))
#' 
#' init.paramt=c(b0t, b1t, log(ht))
#' nll = nll.countST1cov(param=init.paramt, R=Rsites, Tmax=Tsearch,dat=cntm, covar=x)
#' nll
#' 
#' # optimize
#' fit.countST1cov=optim(init.paramt,nll.countST1cov,R=Rsites,Tmax=Tsearch,dat=cntm,covar=x)
#' estpar.countST1cov=fit.countST1cov$par
#' # compare estimates and true parameters
#' c(estpar.countST1cov[1:2],exp(estpar.countST1cov[3]))
#' paramt 
#' 
#' @export
nll.countST1cov=function(param,R,Tmax,dat,covar)
{ 
  ytmat=dat
  b0=param[1]
  b1=param[2]
  h=exp(param[3])
  lam=exp(b0+b1*covar)
  loglik=0
  for(i in 1:R)
  {
    y=ytmat[i,1]
    t1=ytmat[i,2]
    if(y == 0)
    {loglik=loglik-lam[i]*(1-exp(-h*Tmax))}
    if(y >0)
    {
      term1=log(h)+y*log(lam[i])
      term2=(y-1)*log(exp(-h*t1)-exp(-h*Tmax))
      term3=-lam[i]*(1-exp(-h*Tmax))-h*t1
      term4=-log(factorial(y-1))
      loglik=loglik+term1+term2+term3+term4
    }
  }
  return(-loglik)
}