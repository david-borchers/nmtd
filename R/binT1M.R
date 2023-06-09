#' @title Generates time to first detection data for model BinaryT1:M with lambda constant
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site. Then generates an exponential random variable for 
#' each site surveyed, with hazard \eqn{hn}{hn},  where \eqn{h}{h} is the hazard of 
#' detection for a single animal. 
#'
#' @param param A vector comprised of the Poisson rate lambda, and the detection hazard, h.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' 
#' @return Returns an \code{R} by \code{J} matrix of times to first detection values, with 
#' a time equal to \code{Tmax} indicating no detection.
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
#' t1m=as.matrix(generate.binT1M(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' str(t1m)
#' head(t1m)
#' 
#' @export
generate.binT1M=function(param, R, J, Tmax)
{
  lam=param[1]
  h=param[2]
  tmat=matrix(0,R,J)
  for(i in 1:R)
  {
    n=rpois(1,lam)
    if(n==0) {tmat[i,]=Tmax}
    if(n >0) {tvec=rexp(J,h*n)
    tvec[tvec>Tmax]<-Tmax
    tmat[i,]=tvec
    }
  }
  return(tmat)
}

#' @title Evaluates the negative log-likelihood for model BinaryT1:M with lambda constant
#'
#' @description
#' Evaluates the negative log-likelihood for model BinaryT1:M, assuming constant lambda, 
#' given initial parameter estimates and time to first detection data from a 
#' multiple-occasion survey.
#'
#' @param param A vector comprised of the log of the Poisson rate lambda, and the 
#' log of the detection hazard, h.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by \code{J} matrix of time to first detection values, with a 
#' time of \code{Tmax} representing no detection
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
#' t1m=as.matrix(generate.binT1M(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' # optimize
#' init.paramt=c(log(lamt),log(ht))
#' fit.binT1M=optim(init.paramt,nll.binT1M,R=Rsites,J=Jsites,Tmax=Tsearch,dat=t1m)
#' estpar.binT1M=fit.binT1M$par
#' # compare estimates and true parameters
#' exp(estpar.binT1M)
#' paramt 
#' 
#' @export
nll.binT1M=function(param, R, J, Tmax, dat)
{
  lam=exp(param[1])
  h=exp(param[2])
  Jt=matrix(0,R,1)
  W=matrix(0,R,1)
  for(i in 1:R)
  {
    Jt[i]=length(subset(dat[i,], dat[i,] < Tmax))
    # note that W[i] is a sum of all times at site i
    # including those recorded as Tmax, that is ti1 >Tmax
    W[i]=sum(dat[i,])}
  mu=lam*exp(-h*W)
  fin=matrix(0,R,1)
  for(i in 1:R)
  {
    if(Jt[i] > 0)
    {
      term=0
      for(k in 1:Jt[i])
      {
        term=term+(mu[i]^k)*Stirling2(Jt[i],k)}
      fin[i]=fin[i]+log(term)
    }
  }
  loglik=-R*lam+sum(mu)+sum(Jt)*log(h)+
    sum(log(choose(J,Jt)))+sum(fin)
  return(-loglik)
}


#' @title Generates time to first detection data for model Binary:T1M with lambda depending 
#' on a covariate.
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site and depending on the covariate value attached to the site. 
#' Then generates an exponential random variables for each of \code{J} occasions that each
#' site is surveyed, with hazard \eqn{hn}{hn},  where \eqn{h}{h} is the hazard of 
#' detection for a single animal.
#'
#' @param param A vector comprised of parameters \eqn{b_0}{b0}, \eqn{b_1}{b1}, 
#' and the log of the detection hazard (in that order), where the log of the Poisson rate 
#' lambda is equal to \eqn{b_0+b_1x}{b0+b1*x} and \eqn{x}{x} is the covariate \code{covar}.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param covar A vector covariate of length \code{R} on which the expected number of 
#' animals in the site depends linearly (assumed to be the same for all occasions).
#' 
#' @return Returns an \code{R} by \code{J} matrix of times to first detection values, with 
#' a time equal to \code{Tmax} indicating no detection.
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
#' t1m=as.matrix(generate.binT1Mcov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' str(t1m)
#' head(t1m)
#' 
#' @export
generate.binT1Mcov=function(param, R, J, Tmax, covar)
{
  b0=param[1]
  b1=param[2]
  h=param[3]
  lam=exp(b0+b1*covar)
  tmat=matrix(0,R,J)
  for(i in 1:R)
  {
    N=rpois(1,lam[i])
    if(N==0) {tmat[i,]=Tmax}
    if(N >0) {tvec=rexp(J,h*N)
    tvec[tvec>Tmax]<-Tmax
    tmat[i,]=tvec
    }
  }
  return(tmat)
}

#' @title Evaluates the negative log-likelihood for model BinaryT1:M with lambda depending 
#' on a covariate.
#'
#' @description
#' Evaluates the negative log-likelihood for model BinaryT1:M, assuming that lambda 
#' depends on a covariate \code{covar} that is site-dependent, 
#' given initial parameter estimates and binary data from a multiple-occasion survey.
#'
#' @param param A vector comprised of parameters \eqn{b_0}{b0}, \eqn{b_1}{b1}, 
#' and the log of the detection hazard (in that order), where the log of the Poisson rate 
#' lambda is equal to \eqn{b_0+b_1x}{b0+b1*x} and \eqn{x}{x} is the covariate \code{covar}.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by \code{J} matrix of binary values, with a 1 representing
#' detection and a 0 representing no detection
#' @param covar A vector covariate of length \code{R} on which the expected number of 
#' animals in the site depends linearly (assumed to be the same for all occasions).
#' 
#' @return Returns the negative log-likelihood function for model BinaryT1:M evaluated at 
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
#' # data for BinaryT1:M
#' t1m=as.matrix(generate.binT1Mcov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' 
#' init.paramt=c(b0t, b1t, log(ht))
#' nll = nll.binT1Mcov(param=init.paramt, R=Rsites, J=Jsites, Tmax=Tsearch,dat=t1m, covar=x)
#' nll
#' 
#' # optimize
#' fit.binT1Mcov=optim(init.paramt,nll.binT1Mcov,R=Rsites,J=Jsites,Tmax=Tsearch,dat=t1m,covar=x)
#' estpar.binT1Mcov=fit.binT1Mcov$par
#' # compare estimates and true parameters
#' c(estpar.binT1Mcov[1:2],exp(estpar.binT1Mcov[3]))
#' paramt 
#' 
#' @export
nll.binT1Mcov=function(param, R, J, Tmax, dat, covar)
{
  b0=param[1]
  b1=param[2]
  h=exp(param[3])
  lam=exp(b0+b1*covar)
  Jt=matrix(0,R,1)
  W=matrix(0,R,1)
  for(i in 1:R)
  {
    Jt[i]=length(subset(dat[i,], dat[i,] < Tmax))
    # note that W[i] is a sum of all times at site i
    # including those recorded as Tmax, that is ti1 >Tmax
    W[i]=sum(dat[i,])
  }
  mu=lam*exp(-h*W)
  fin=matrix(0,R,1)
  for(i in 1:R)
  {
    if(Jt[i] > 0)
    {
      term=0
      for(k in 1:Jt[i])
      {
        term=term+(mu[i]^k)*Stirling2(Jt[i],k)}
      fin[i]=fin[i]+log(term)
    }
  }
  # note that the term -lam_i (1-exp(-h W_i)) = -lam_i+[lam_i exp(-h W_i)]
  # covers all R since if Jt[i]=0 then W[i] = h*Tmax*J
  loglik=-sum(lam)+sum(mu)+sum(Jt)*log(h)+ sum(log(choose(J,Jt)))+sum(fin)
  return(-loglik)
}
