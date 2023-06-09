#' @title Generates count data for model CountT:S with lambda constant
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site. 
#' Then generates binomial random variables for each 
#' site on each of \code{J} occasions, given the number of animals present, 
#' assuming a constant hazard of detection for each animal. 
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
#' Jsites=5    #multiple visits
#' Tsearch=3 #maximum time

#' #true parameters
#' lamt = 2
#' ht=0.4620981
#' paramt=c(lamt,ht)
#' 
#' # data for Count:M
#' cntm=as.matrix(generate.countM(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' str(cntm)
#' head(cntm)
#' 
#' @export
generate.countM=function(param, R, J, Tmax)
{
  lam=param[1]
  h=param[2]
  ymat=matrix(0,R,J)
  for(i in 1:R){
    n=rpois(1,lam)
    p=1-exp(-h*Tmax)
    for(j in 1:J)
    {
      ymat[i,j]=rbinom(1,n,p)}
  }
  return(ymat)
}

#' @title Evaluates the negative log-likelihood for model Count:M with lambda constant
#'
#' @description
#' Evaluates the negative log-likelihood for model Count:M, assuming constant lambda, 
#' given initial parameter estimates and count data from a multiple-occasion survey.
#'
#' @param param A vector comprised of the log of the Poisson rate lambda, and the 
#' log of the detection hazard, h.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by \code{J} matrix of counts.
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
#' # data for Count:M
#' set.seed(123) # for reproducibility
#' cntm=as.matrix(generate.countM(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' # optimize
#' init.paramt=c(log(lamt),log(ht))
#' fit.countM=optim(init.paramt,nll.countM,R=Rsites,J=Jsites,Tmax=Tsearch,dat=cntm)
#' estpar.countM=fit.countM$par
#' # compare estimates and true parameters
#' exp(estpar.countM)
#' paramt 
#' 
#' @export
nll.countM=function(param, R, J, Tmax, dat)
{
  lam=exp(param[1])
  h=exp(param[2])
  loglik=0
  p=1-exp(-h*Tmax)
  ymat=dat
  for(i in 1:R)
  {
    # set up yvec, p, lograt
    yvec=ymat[i,]
    lograt=sum(yvec)*log(p/(1-p))
    #sort yvec
    yvec=sort(yvec)
    yij=yvec[-J]
    yiJ=yvec[J]
    term=sum(log(choose(yiJ,yij)))
    # calculate hypergeometric function
    theta=lam*(1-p)^J
    hgfun=genhypergeo(rep(yiJ,(J-1)) + 1,yiJ - yij + 1,theta)
    loglik=loglik-lam+lograt+yiJ*log(theta)+
      log(hgfun)+term-log(factorial(yiJ))
  }
  return(-loglik)
}


#' @title Generates binary data for model Count:M with lambda depending 
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
#' 1-exp(-ht*Tmax)      #P(individual detected)
#' 1-exp(-mean(lambda)) #P(site occupancy)
#' paramt=c(b0t, b1t,ht)
#' 
#' # data for Count:M
#' cntm=as.matrix(generate.countMcov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' str(cntm)
#' head(cntm)
#' 
#' @export
generate.countMcov=function(param,R,J,Tmax,covar)
{
  b0=param[1]
  b1=param[2]
  lam=exp(b0+b1*covar)
  h=param[3]
  ymat=matrix(0,R,J)
  for(i in 1:R){
    n=rpois(1,lam[i])
    p=1-exp(-h*Tmax)
    for(j in 1:J)
    {
      ymat[i,j]=rbinom(1,n,p)}
  }
  return(ymat)
}

#' @title Evaluates the negative log-likelihood for model Count:M with lambda depending 
#' on a covariate.
#'
#' @description
#' Evaluates the negative log-likelihood for model Count:M, assuming that lambda 
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
#' x=rnorm(Rsites,0,1)
#' lambda=exp(b0t+b1t*x)
#' mean(lambda)
#' 
#' # true probabilites
#' 1-exp(-ht*Tmax)      #P(individual detected)
#' 1-exp(-mean(lambda)) #P(site occupancy)
#' paramt=c(b0t,b1t,ht)
#' 
#' # data for Count:M
#' cntm=as.matrix(generate.countMcov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' 
#' init.paramt=c(b0t, b1t, log(ht))
#' nll = nll.countMcov(param=init.paramt, R=Rsites, J=Jsites, Tmax=Tsearch,dat=cntm, covar=x)
#' nll
#' 
#' # optimize
#' fit.countMcov=optim(init.paramt,nll.countMcov,R=Rsites,J=Jsites,Tmax=Tsearch,dat=cntm,covar=x)
#' estpar.countMcov=fit.countMcov$par
#' # compare estimates and true parameters
#' c(estpar.countMcov[1:2],exp(estpar.countMcov[3]))
#' paramt 
#' 
#' @export
nll.countMcov=function(param,R,J,Tmax,dat,covar)
{
    b0=param[1]
    b1=param[2]
    h=exp(param[3])
    lam=exp(b0+b1*covar)
    loglik=0
    ymat=dat
    for(i in 1:R)
    {
      # set up yvec, p, lograt, theta
      yvec=ymat[i,]
      p=1-exp(-h*Tmax)
      lograt=sum(yvec*log(p/(1-p)))
      theta=lam[i]*(1-p)^J
      # now sort yvec
      yvec=sort(yvec)
      yij=yvec[-J]
      yiJ=yvec[J]
      term=sum(log(choose(yiJ,yij)))
      hgfun=genhypergeo(rep(yiJ,(J-1)) + 1,yiJ - yij + 1,theta)
      # now calculate log-likelihood
      loglik=loglik+lograt-lam[i]+yiJ*log(theta)-
        log(factorial(yiJ))+term+log(hgfun)
    }
    return(-loglik)
  }