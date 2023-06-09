#' @title Generates count data for model CountT:M with lambda constant
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site on each of \code{J} occasions. 
#' Then generates exponential random variables for each animal in each site, 
#' assuming a constant hazard of detection for each animal. 
#'
#' @param param A vector comprised of the Poisson rate lambda, and the detection hazard, h.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' 
#' @return Returns an \code{R} by \code{2J} matrix with first \code{J} columns being the numbers 
#' of detections in each site on each occasion and the remaining \code{J} columns being the 
#' **sum** of times to detection for all detected  animals within the site on the occasion 
#' (with zero representing no detections).
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
#' # data for CountT:M
#' cntm=as.matrix(generate.countMT(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' str(cntm)
#' head(cntm)
#' 
#' @export
generate.countMT=function(param, R, J, Tmax)
{
  lam=param[1]
  h=param[2]
  ymat=matrix(0,R,J)
  tsmat=matrix(0,R,J)
  for(i in 1:R)
  {
    n=rpois(1,lam)
    if(n==0){ymat[i,]=matrix(0,J,1)
    tsmat[i,]=matrix(0,J,1)}
    if( n > 0)
      for(j in 1:J)
      {
        tvec=sort(rexp(n,h))
        tvecltT=as.matrix(tvec[tvec < Tmax])
        ymat[i,j]=nrow(tvecltT)
        tsmat[i,j]=sum(tvecltT)
      }
  }
  return(cbind(ymat,tsmat))
}

#' @title Evaluates the negative log-likelihood for model CountT:M with lambda constant
#'
#' @description
#' Evaluates the negative log-likelihood for model CountT:M, assuming constant lambda, 
#' given initial parameter estimates and count data from a multiple-occasion survey.
#'
#' @param param A vector comprised of the log of the Poisson rate lambda, and the 
#' log of the detection hazard, h.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by \code{2J} matrix with first \code{J} columns being the numbers 
#' of detections in each site on each occasion and the remaining \code{J} columns being the 
#' **sum** of times to detection for all detected  animals within the site on the occasion 
#' (with zero representing no detections).
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
#' # data for CountT:M
#' set.seed(123) # for reproducibility
#' cntm=as.matrix(generate.countMT(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' # optimize
#' init.paramt=c(log(lamt),log(ht))
#' fit.countMT=optim(init.paramt,nll.countMT,R=Rsites,J=Jsites,Tmax=Tsearch,dat=cntm)
#' estpar.countMT=fit.countMT$par
#' # compare estimates and true parameters
#' exp(estpar.countMT)
#' paramt 
#' 
#' @export
nll.countMT=function(param, R, J, Tmax, dat)
{ 
  ymat=dat[,1:J]
  tsmat=dat[,(J+1):(2*J)]
  lam=exp(param[1])
  h=exp(param[2])
  loglik=0
  for(i in 1:R)
  {
    yvec=ymat[i,]
    yvec=sort(yvec, decreasing=FALSE)
    yij=yvec[-J]
    yiJ=yvec[J]
    tsum=tsmat[i,]
    term1=sum(yvec)*log(h)-h*sum(tsum)
    term2=-lam-h*sum(Tmax*(yiJ-yvec))+yiJ*log(lam)
    term3=(J-1)*log(factorial(yiJ))-sum(log(factorial(yiJ-yij)))
    theta=lam*exp(-h*Tmax*J)
    hgfun=genhypergeo(rep(yiJ,(J-1)) + 1,yiJ - yij + 1,theta)
    term4=log(hgfun)
    loglik=loglik+term1+term2+term3+term4
  }
  return(-loglik)
}

#' @title Generates binary data for model CountT:M with lambda depending 
#' on a covariate.
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site and depending on the covariate value attached to the site. 
#' Then generates Poissin random variables for each animal on each of \code{J} occasions that 
#' each site is surveyed, assuming a constant hazard of detection for each animal. 
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
#' @return Returns an \code{R} by \code{2J} matrix with first \code{J} columns being the numbers 
#' of detections in each site on each occasion and the remaining \code{J} columns being the 
#' **sum** of times to detection for all detected  animals within the site on the occasion 
#' (with zero representing no detections).
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
#' # data for CountT:M
#' cntm=as.matrix(generate.countMTcov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' str(cntm)
#' head(cntm)
#' 
#' @export
generate.countMTcov=function(param,R,J,Tmax,covar)
{
  b0=param[1]
  b1=param[2]
  h=param[3]
  lam=exp(b0+b1*covar)
  ymat=matrix(0,R,J)
  tsmat=matrix(0,R,J)
  for(i in 1:R)
  { n=rpois(1,lam[i])
  if( n > 0)
    for(j in 1:J)
    {
      tvec=rexp(n,h)
      tvecltT=as.matrix(tvec[tvec < Tmax])
      ymat[i,j]=nrow(tvecltT)
      tsmat[i,j]=sum(tvecltT)
    }
  }
  return(cbind(ymat,tsmat))
}


#' @title Evaluates the negative log-likelihood for model CountT:M with lambda depending 
#' on a covariate.
#'
#' @description
#' Evaluates the negative log-likelihood for model CountT:M, assuming that lambda 
#' depends on a covariate \code{covar} that is site-dependent, 
#' given initial parameter estimates and binary data from a multiple-occasion survey.
#'
#' @param param A vector comprised of parameters \eqn{b_0}{b0}, \eqn{b_1}{b1}, 
#' and the log of the detection hazard (in that order), where the log of the Poisson rate 
#' lambda is equal to \eqn{b_0+b_1x}{b0+b1*x} and \eqn{x}{x} is the covariate \code{covar}.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by \code{2J} matrix with first \code{J} columns being the numbers 
#' of detections in each site on each occasion and the remaining \code{J} columns being the 
#' **sum** of times to detection for all detected  animals within the site on the occasion 
#' (with zero representing no detections).
#' @param covar A vector covariate of length \code{R} on which the expected number of 
#' animals in the site depends linearly (assumed to be the same for all occasions).
#' 
#' @return Returns the negative log-likelihood function for model CountT:M evaluated at 
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
#' # data for CountT:M
#' cntm=as.matrix(generate.countMTcov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' 
#' init.paramt=c(b0t, b1t, log(ht))
#' nll = nll.countMTcov(param=init.paramt, R=Rsites, J=Jsites, Tmax=Tsearch,dat=cntm, covar=x)
#' nll
#' 
#' # optimize
#' fit.countMTcov=optim(init.paramt,nll.countMTcov,R=Rsites,J=Jsites,Tmax=Tsearch,dat=cntm,covar=x)
#' estpar.countMTcov=fit.countMTcov$par
#' # compare estimates and true parameters
#' c(estpar.countMTcov[1:2],exp(estpar.countMTcov[3]))
#' paramt 
#' 
#' @export
nll.countMTcov=function(param,R,J,Tmax,dat,covar)
{ 
  ymat=dat[,1:J]
  tsmat=dat[,(J+1):(2*J)]
  b0=param[1]
  b1=param[2]
  h=exp(param[3])
  lam=exp(b0+b1*covar)
  loglik=0
  for(i in 1:R)
  {
    yvec=ymat[i,]
    tsum=tsmat[i,]
    yvec=sort(yvec, decreasing=FALSE)
    yij=yvec[-J]
    yiJ=yvec[J]
    term1=sum(yvec)*log(h)-h*sum(tsum)
    term2=-lam[i]-h*sum(Tmax*(yiJ-yvec))+yiJ*log(lam[i])
    term3=(J-1)*log(factorial(yiJ))-sum(log(factorial(yiJ-yij)))
    theta=lam[i]*exp(-h*Tmax*J)
    hgfun=genhypergeo(rep(yiJ,(J-1)) + 1,yiJ - yij + 1,theta)
    term4=log(hgfun)
    loglik=loglik+term1+term2+term3+term4}
  return(-loglik)
}
