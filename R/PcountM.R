#' @title Generates count data for model PCount:M with lambda constant
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site. Then generates Poisson random variables for each 
#' site, on each of \code{J} occassion, given the number of animals present, assuming a 
#' constant hazard of detection for each animal. 
#'
#' @param param A vector comprised of the Poisson rate lambda, and the detection hazard, h.
#' @param R The number of sites.
#' @param J The number of occasions.
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' 
#' @return Returns an \code{R} by \code{J} matrix of counts.
#' 
#' @examples 
#' # seed
#' set.seed(123)
#' # setting
#' Rsites=100. #number of sites
#' Tsearch=3 #maximum time
#' Joccs=5 #number of occasions
#' #true parameters
#' lamt=2
#' gt=0.4620981
#' paramt=c(lamt,gt)
#' 
#' ymat=generate.PcountM(paramt, R=Rsites, J=Joccs, Tmax=Tsearch)
#' str(ymat)
#' head(ymat)
#' 
#' @export
generate.PcountM = function(param, R, J, Tmax)
{
  lam=param[1]
  g=param[2]
  ymat=matrix(0,R,J)
  for(i in 1:R)
  {
    n=rpois(1,lam)
    for(j in 1:J) ymat[i,j]=rpois(1,g*n*Tmax)
  }
  return(ymat)
}


#' @title Evaluates the negeative log-likelihood for model PCount:M with lambda constant.
#'
#' @description
#' Evaluates the negeative log-likelihood for model PCount:M, assuming that lambda 
#' is constant, given initial parameter estimates and binary data from a multiple-occasion 
#' survey.
#'
#' @param param A vector comprised of the Poisson rate lambda, and the detection hazard, h.
#' @param R The number of sites.
#' @param J The number of occasions.
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by \code{J} matrix of counts
#' 
#' @return Returns the negative log-likelihood function evaluated for model PCount:M at 
#' the parameter values passed in \code{param}.
#' 
#' @examples 
#' # seed
#' set.seed(123)
#' # setting
#' Rsites=100. #number of sites
#' Tsearch=3 #maximum time
#' Joccs=5 #number of occasions
#' #true parameters
#' lamt=2
#' gt=0.4620981
#' paramt=c(lamt,gt)
#' 
#' ymat=generate.PcountM(paramt, R=Rsites, J=Joccs, Tmax=Tsearch)
#' str(ymat)
#' head(ymat)
#' 
# starting values for optimiser
#' param0=c(log(lamt),log(gt))
#' # evaluate negative log-likelihood at starting values:
#' nll.PcountM(param0, R=Rsites, J=Joccs, Tmax=Tsearch, dat=ymat)
#' 
#' # optimize
#' fit.PcountM=optim(param0, nll.PcountS, R=Rsites, J=Joccs, Tmax=Tsearch, 
#'                   dat=ymat, method="BFGS")
#' estpar.PcountM=fit.PcountM$par
#' # compare estimates and true parameters
#' exp(estpar.PcountM)
#' paramt 
#' 
#' @export
nll.PcountM = function(param, R, J, Tmax, dat)
{
  lam=exp(param[1])
  g=exp(param[2])
  mu=lam*exp(-g*Tmax*J)
  fin=matrix(0,R,1)
  term=matrix(0,R,1)
  ytot=rowSums(ymat)
  for(i in 1:R)
  { 
    if(ytot[i] > 0)
    {
      for(k in 0:ytot[i])
        term[i]=term[i]+exp(k*log(mu)+log(Stirling2(ytot[i],k)))
      fin[i]=fin[i]+log(term[i])
    }
  }
  loglik=-R*lam+R*mu+sum(ymat)*log(g*Tmax)+sum(fin)
  return(-loglik)
}


#' @title Generates count data for model PCount:M with lambda depending on a covariate.
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site and depending on the covariate value attached to the site. 
#' Then generates Poisson random variables for each site, oneach of \code{J} occasions, 
#' given the number of animals present, assuming a constant hazard of detection for each animal. 
#'
#' @param param A vector comprised of the Poisson rate lambda, and the detection hazard, h.
#' @param R The number of sites.
#' @param J The number of occasions
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param covar A vector covariate of length \code{R} on which the expected number of 
#' animals in the site depends linearly (assumed to be the same for all occasions).
#' 
#' @return Returns an \code{R} by \code{J} matrix of counts.
#' 
#' @examples
#' # seed
#' set.seed(321)
#' # setting
#' Rsites=100. #number of sites
#' Joccs=5 # number of occasions
#' Tsearch=3 #maximum time
#' #true parameters
#' b0t =0.660878  
#' b1t =0.255413
#' gt=0.4620981
#' # simulation of the covariate x over R sites for true parameters
#' x=rnorm(R,0,1)
#' hist(x)
#' lambda=exp(b0t+b1t*x)
#' mean(lambda)
#' # true parameters
#' paramt=c(b0t,b1t,gt)
#' 
#' # generate
#' ymat=generate.PcountMcov(paramt, R=Rsites, J=Joccs, Tmax=Tsearch, covar=x)
#' str(ymat)
#' head(ymat)
#' 
#' @export
generate.PcountMcov = function(param, R, J, Tmax, covar)
{
  b0=param[1]
  b1=param[2]
  g=param[3]
  lam=exp(b0+b1*covar)
  ymat=matrix(0,R,J)
  for(i in 1:R)
  {n=rpois(1,lam[i])
  for(j in 1:J)
    ymat[i,j]=rpois(1,g*n*Tmax)
  }
  return(ymat)
}


#' @title Evaluates the negeative log-likelihood for model PCount:M with lambda depending 
#' on a covariate.
#'
#' @description
#' Evaluates the negeative log-likelihood for model PCount:M, assuming that lambda 
#' depends on a covariate \code{covar} that is site-dependent, 
#' given initial parameter estimates and Poisson count data from a \code{J}-occasion survey.
#'
#' @param param A vector comprised of parameters \eqn{b_0}{b0}, \eqn{b_1}{b1}, 
#' and the log of the detection hazard (in that order), where the log of the Poisson rate 
#' lambda is equal to \eqn{b_0+b_1x}{b0+b1*x} and \eqn{x}{x} is the covariate \code{covar}.
#' @param R The number of sites.
#' @param J The number of occasions
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by 1 matrix of counts
#' @param covar A vector covariate of length \code{R} on which the expected number of 
#' animals in the site depends linearly (assumed to be the same for all occasions).
#' 
#' @return Returns the negative log-likelihood function for model PCount:M evaluated at 
#' the parameter values passed in \code{param}.
#' 
#' @examples
#' # seed
#' set.seed(321)
#' # setting
#' Rsites=100. #number of sites
#' Joccs=5 # number of occasions
#' Tsearch=3 #maximum time
#' #true parameters
#' b0t =0.660878  
#' b1t =0.255413
#' gt=0.4620981
#' # simulation of the covariate x over R sites for true parameters
#' x=rnorm(R,0,1)
#' lambda=exp(b0t+b1t*x)
#' mean(lambda)
#' # true parameters
#' paramt=c(b0t,b1t,gt)
#' 
#' # generate
#' ymat=generate.PcountMcov(paramt, R=Rsites, J=Joccs, Tmax=Tsearch, covar=x)
#' str(ymat)
#' head(ymat)
#' 
#' # starting values 
#' param0=cbind(b0t,b1t,log(gt))
#' # evaluate likelihood at starting value
#' nll.PcountMcov(param0, R=Rsites, J=Joccs, Tmax=Tsearch, dat=ymat, covar=x)
#' 
#' # optimize
#' fit.PcountMcov=optim(param0, nll.PcountMcov, R=Rsites, J=Joccs, Tmax=Tsearch, 
#'                      dat=ymat,covar=x,,method="BFGS")
#' estpar.PcountMcov=fit.PcountMcov$par
#' # compare estimates and true parameters
#' c(estpar.PcountMcov[1:2],exp(estpar.PcountMcov[3]))
#' paramt 
#' 
#' @export
nll.PcountMcov = function(param, R, J, Tmax, dat, covar)
{
  b0=param[1]
  b1=param[2]
  g=exp(param[3])
  lam=exp(b0+b1*covar)
  fin=matrix(0,R,1)
  term=matrix(0,R,1)
  ytot=rowSums(ymat)
  mu=lam*exp(-g*Tmax*J)
  for(i in 1:R)
  {
    if(ytot[i] > 0)
    {
      for(k in 0:ytot[i])
        term[i]=term[i]+exp(k*log(mu[i])+log(Stirling2(ytot[i],k)))
      fin[i]=fin[i]+log(term[i])}
  }
  loglik=-sum(lam)+sum(mu)+sum(ymat)*log(g*Tmax)+sum(fin)
  return(-loglik)
}
