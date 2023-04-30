#' @title Generates binary data for model BinaryT1:S with constant lambda.
#'
#' @description
#' Generates Poisson random variables \eqn{n}{n} for each of the \code{R} sites, corresponding to the 
#' number of animals in each site. 
#' Then generates an exponential random variable for each site surveyed, with hazard \eqn{hn}{hn}, 
#' where \eqn{h}{h} is the hazard of detection for a single animal. 
#'
#' @param param A vector comprised of the Poisson rate parameter lambda, and the detection hazard
#' \eqn{h}{h}.
#' @param R The number of sites.
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' 
#' @return Returns an \code{R} by 1 matrix containing the times of
#' first detections in each site, with a time equal to \code{Tmax} indication no detection.
#' 
#' @examples 
#' Rsites=100. #number of sites
#' Tsearch=3 #maximum time
#' 
#' #true parameters
#' lamt = 2
#' ht=0.4620981
#' paramt=c(lamt,ht)
#' 
#' # data for BinaryT1:S
#' bns=as.matrix(generate.binT1S(paramt,R=Rsites,Tmax=Tsearch))
#' str(bns)
#' head(bns)
#' 
#' @export
generate.binT1S=function(param,R,Tmax)
{
  lam=param[1]
  h=param[2]
  tvec=matrix(0,R,1)
  for(i in 1:R)
  {
    n=rpois(1,lam);
    if(n==0)
    {tvec[i]=Tmax}
    if(n > 0)
    {ti=rexp(1,h*n);
    if(ti>Tmax)
    {tvec[i]=Tmax}
    if(ti<Tmax)
    {tvec[i]=ti}
    }
  }
  return(tvec)
}

#' @title Evaluates the negeative log-likelihood for model BinaryT1:S with lambda depending 
#' on a covariate.
#'
#' @description
#' Evaluates the negeative log-likelihood for model BinaryT1:S, assuming that lambda 
#' is constant, given initial parameter estimates and binary data from a 
#' single-occasion survey.
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
#' @return Returns the negative log-likelihood function evaluated for model BinaryT1:S at 
#' the parameter values passed in \code{param}.
#' 
#' @examples 
#' Rsites=100. #number of sites
#' Tsearch=3 #maximum time
#' 
#' #true parameters
#' #true parameters
#' lamt = 2
#' ht=0.4620981
#' paramt=c(lamt,ht)
#' 
#' # data for BinaryT1:S
#' bns=as.matrix(generate.binT1S(paramt,R=Rsites,Tmax=Tsearch))
#' 
#' init.paramt=c(log(lamt),log(ht))
#' nll = nll.binT1S(param=init.paramt, R=Rsites, Tmax=Tsearch,dat=bns)
#' nll
#' 
#' # optimize
#' set.seed(123) # for reproducibility
#' bn1s=as.matrix(generate.binT1S(paramt,R=Rsites,Tmax=Tsearch))
#' # optimize
#' fit.bin1S=optim(init.paramt,nll.binT1S,R=Rsites,Tmax=Tsearch,dat=bn1s)
#' estpar.bin1S=fit.bin1S$par
#' # compare estimates and true parameters
#' exp(estpar.bin1S)
#' paramt 
#' 
#' @export
nll.binT1S=function(param,R,Tmax,dat)
{
  lam=exp(param[1])
  h=exp(param[2])
  ti=as.matrix(subset(dat,dat < Tmax))
  Rt=nrow(ti) # total number of sites at which t_i1 < Tmax
  RT=R-Rt
  sumt=sum(ti)
  loglik=-RT*lam*(1-exp(-h*Tmax))+
    Rt*log(h*lam) -lam*sum(1-exp(-h*ti))-h*sumt
  return(-loglik)
}

#' @title Generates binary data for model BinaryT1:S with lambda depending 
#' on a covariate.
#'
#' @description
#' Generates Poisson random variables \eqn{n}{n} for each of the \code{R} sites, corresponding to the 
#' number of animals in each site. 
#' Then generates an exponential random variable for each site surveyed, with hazard \eqn{hn}{hn}, 
#' where \eqn{h}{h} is the hazard of detection for a single animal. 
#'
#' @param param A vector comprised of parameters \eqn{b_0}{b0}, \eqn{b_1}{b1}, 
#' and the log of the detection hazard \eqn{h}{h} (in that order), where the log of the Poisson rate 
#' lambda is equal to \eqn{b_0+b_1x}{b0+b1*x} and \eqn{x}{x} is the covariate \code{covar}.
#' @param R The number of sites.
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param covar A vector covariate of length \code{R} on which the expected number of 
#' animals in the site depends linearly (assumed to be the same for all occasions).
#' 
#' @return Returns an \code{R} by 1 matrix containing the times of first detections in 
#' each site, with a time equal to \code{Tmax} indication no detection.
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
#' x=rnorm(R,0,1)
#' lambda=exp(b0t+b1t*x)
#' mean(lambda)
#' 
#' # true probabilites
#' 1-exp(-ht*Tmax)      #P(individual detected)
#' 1-exp(-mean(lambda)) #P(site occupancy)
#' paramt=c(b0t, b1t,ht)
#' 
#' # data for BinaryT1:Scov
#' bns=as.matrix(generate.binT1Scov(paramt,R=Rsites,Tmax=Tsearch,covar=x))
#' str(bns)
#' head(bns)
#' 
#' @export
generate.binT1Scov=function(param,R,Tmax,covar)
{
  b0=param[1]
  b1=param[2]
  h=param[3]
  lam=exp(b0+b1*covar)
  tvec=matrix(0,R,1)
  for(i in 1:R)
  {
    n=rpois(1,lam[i]);
    if(n==0)
    {tvec[i]=Tmax}
    if(n > 0)
    {ti=rexp(1,h*n);
    if(ti>Tmax)
    {tvec[i]=Tmax}
    if(ti<Tmax)
    {tvec[i]=ti}
    }
  }
  return(tvec)
}

#' @title Evaluates the negeative log-likelihood for model BinaryT1:S with lambda depending 
#' on a covariate.
#'
#' @description
#' Evaluates the negeative log-likelihood for model BinaryT1:S, assuming that lambda 
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
#' @return Returns the negative log-likelihood function evaluated for model BinaryT1:S at 
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
#' x=rnorm(R,0,1)
#' lambda=exp(b0t+b1t*x)
#' mean(lambda)
#' 
#' # true probabilites
#' 1-exp(-ht*Tmax)      #P(individual detected)
#' 1-exp(-mean(lambda)) #P(site occupancy)
#' paramt=c(b0t,b1t,ht)
#' 
#' # data for BinaryT1:Scov
#' bns=as.matrix(generate.binT1Scov(paramt,R=Rsites,Tmax=Tsearch,covar=x))
#' 
#' init.paramt=c(b0t, b1t, log(ht))
#' nll = nll.binT1Scov(param=init.paramt, R=Rsites, Tmax=Tsearch,dat=bns, covar=x)
#' nll
#' 
#' # optimize
#' fit.binT1Scov=optim(init.paramt,nll.binT1Scov,R=Rsites,Tmax=Tsearch,dat=bns,covar=x)
#' estpar.binT1Scov=fit.binT1Scov$par
#' # compare estimates and true parameters
#' c(estpar.binT1Scov[1:2],exp(estpar.binT1Scov[3]))
#' paramt 
#' 
#' @export
nll.binT1Scov=function(param,R,Tmax,dat,covar)
{
  b0=param[1]
  b1=param[2]
  h=exp(param[3])
  lam=exp(b0+b1*covar)
  loglikt=0
  loglikT=0
  for(i in 1:R)
  {
    if(dat[i] < Tmax)
    {loglikt=loglikt+log(h*lam[i])-lam[i]+ lam[i]*exp(-h*dat[i])-
      h*dat[i]}
    if(dat[i] == Tmax)
    {loglikT=loglikT-lam[i]*(1-exp(-h*Tmax))}
  }
  loglik=sum(loglikt+loglikT)
  return(-loglik)
}

