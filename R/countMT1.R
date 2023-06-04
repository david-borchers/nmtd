#' @title Generates count data for model CountT1:M with lambda constant
#'
#' @description
#' Generates Poisson random variables for each of the \code{R} sites, corresponding to the 
#' number of animals in each site on each of \code{J} occasions. 
#' Then generates exponential random variables for each ainimal in each site, 
#' assuming a constant hazard of detection for each animal. 
#'
#' @param param A vector comprised of the Poisson rate lambda, and the detection hazard, h.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' 
#' @return Returns an \code{R} by \code{2J} matrix with first \code{J} columns being the numbers 
#' of detections in each site on each occasion and the remaining \code{J} columns being the 
#' shortest time to detection for all detected  animals with in the site on the occasion 
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
#' # data for CountT1:M
#' cntm=as.matrix(generate.countMT1(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' str(cntm)
#' head(cntm)
#' 
#' @export
generate.countMT1=function(param, R, J, Tmax)
{
  lam=param[1]
  h=param[2]
  ymat=matrix(0,R,J)
  t1mat=matrix(0,R,J)
  for(i in 1:R)
  { n=rpois(1,lam)
  if( n > 0)
    for(j in 1:J)
    {
      tvec=sort(rexp(n,h))
      tvecltT=as.matrix(tvec[tvec < Tmax])
      y=nrow(tvecltT)
      if(y==0)
      {ymat[i,j]=0
      t1mat[i,j]=0}else
      {ymat[i,j]=y
      t1mat[i,j]=min(tvecltT)}
    }
  }
  return(cbind(ymat,t1mat))
}



#' @title Evaluates the negeative log-likelihood for model CountT:M with lambda constant
#'
#' @description
#' Evaluates the negeative log-likelihood for model CountT:M, assuming constant lambda, 
#' given initial parameter estimates and count data from a multiple-occasion survey.
#'
#' @param param A vector comprised of the log of the Poisson rate lambda, and the 
#' log of the detection hazard, h.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by \code{2J} matrix with first \code{J} columns being the numbers 
#' of detections in each site on each occasion and the remaining \code{J} columns being the 
#' times to first detection of an animal within the site on the occasion 
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
#' # data for CountT1:M
#' set.seed(123) # for reproducibility
#' cntm=as.matrix(generate.countMT1(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
#' # optimize
#' init.paramt=c(log(lamt),log(ht))
#' fit.countMT1=optim(init.paramt,nll.countMT1,R=Rsites,J=Jsites,Tmax=Tsearch,dat=cntm)
#' estpar.countMT1=fit.countMT1$par
#' # compare estimates and true parameters
#' exp(estpar.countMT1)
#' paramt 
#' 
#' @export
nll.countMT1=function(param, R, J, Tmax, dat)
{ 
  ymat=dat[,1:J]
  t1mat=dat[,(J+1):(2*J)]
  lam=exp(param[1])
  h=exp(param[2])
  # all R
  loglik=0
  for(i in 1:R)
  {
    yvec=ymat[i,]
    if(sum(yvec)==0)
    {
      loglik=loglik-lam*(1-exp(-h*Tmax*J))
    }
    if(sum(yvec)>0)
    {# sort on yvec and t1
      t1vec=t1mat[i,]
      yt=cbind(yvec,t1vec)
      yt=yt[order(yt[,1],decreasing=FALSE),]
      yvec=yt[,1]
      t1vec=yt[,2]
      yij=yvec[-J]
      ymax=yvec[J]
      # terms for y[i,j] > 0 in site i
      for(j in 1:J)
      {
        y=yvec[j]
        if(y > 0)
        {
          term1=log(h)-h*t1vec[j]
          term2=(y-1)*log((exp(-h*t1vec[j])-exp(-h*Tmax)))
          term3=-sum(log(factorial(y-1)))
          loglik=loglik+term1+term2+term3
        }
      }
      #terms for all J equiv to sum over n
      term1=ymax*log(lam)
      term2=-lam-h*Tmax*sum(ymax-yvec)
      term3=(J-1)*log(factorial(ymax))-sum(log(factorial(ymax-yij)))
      theta=lam*exp(-h*Tmax*J)
      hgfun=genhypergeo(rep(ymax,(J-1)) + 1,ymax - yij + 1,theta)
      term4=log(hgfun)
      loglik=loglik+term1+term2+term3+term4}
  }
  return(-loglik)
}




#' @title Generates binary data for model CountT1:M with lambda depending 
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
#' times to first detection of an animal within the site on the occasion 
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
#' cntm=as.matrix(generate.countMT1cov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' str(cntm)
#' head(cntm)
#' 
#' @export
generate.countMT1cov=function(param,R,J,Tmax,covar)
{
  b0=param[1]
  b1=param[2]
  h=param[3]
  lam=exp(b0+b1*covar)
  ymat=matrix(0,R,J)
  t1mat=matrix(0,R,J)
  for(i in 1:R)
  { n=rpois(1,lam[i])
  if( n > 0)
    for(j in 1:J)
    {
      tvec=sort(rexp(n,h))
      tvecltT=as.matrix(tvec[tvec < Tmax])
      y=nrow(tvecltT)
      if(y>0)
      {ymat[i,j]=y
      t1mat[i,j]=min(tvecltT)}
    }
  }
  return(cbind(ymat,t1mat))
}




#' @title Evaluates the negeative log-likelihood for model CountT1:M with lambda depending 
#' on a covariate.
#'
#' @description
#' Evaluates the negeative log-likelihood for model CountT1:M, assuming that lambda 
#' depends on a covariate \code{covar} that is site-dependent, 
#' given initial parameter estimates and binary data from a multiple-occasion survey.
#'
#' @param param A vector comprised of parameters \eqn{b_0}{b0}, \eqn{b_1}{b1}, 
#' and the log of the detection hazard (in that order), where the log of the Poisson rate 
#' lambda is equal to \eqn{b_0+b_1x}{b0+b1*x} and \eqn{x}{x} is the coaviate \code{covar}.
#' @param R The number of sites.
#' @param J The number of occasions (assumed the same for all sites).
#' @param Tmax The survey duration (assumed to be the same for all sites)
#' @param dat An \code{R} by \code{2J} matrix with first \code{J} columns being the numbers 
#' of detections in each site on each occasion and the remaining \code{J} columns being the 
#' times to first detection of an animal within the site on the occasion 
#' (with zero representing no detections).
#' @param covar A vector covariate of length \code{R} on which the expected number of 
#' animals in the site depends linearly (assumed to be the same for all occasions).
#' 
#' @return Returns the negative log-likelihood function for model CountT1:M evaluated at 
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
#' cntm=as.matrix(generate.countMT1cov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
#' 
#' init.paramt=c(b0t, b1t, log(ht))
#' nll = nll.countMT1cov(param=init.paramt, R=Rsites, J=Jsites, Tmax=Tsearch,dat=cntm, covar=x)
#' nll
#' 
#' # optimize
#' fit.countMT1cov=optim(init.paramt,nll.countMT1cov,R=Rsites,J=Jsites,Tmax=Tsearch,dat=cntm,covar=x)
#' estpar.countMT1cov=fit.countMT1cov$par
#' # compare estimates and true parameters
#' c(estpar.countMT1cov[1:2],exp(estpar.countMT1cov[3]))
#' paramt 
#' 
#' @export
nll.countMT1cov=function(param,R,J,Tmax,dat,covar)
{
  b0=param[1]
  b1=param[2]
  h=exp(param[3])
  ymat=dat[,1:J]
  t1mat=dat[,(J+1):(2*J)]
  lam=exp(b0+b1*covar)
  loglik=0
  for(i in 1:R)
  {
    yvec=ymat[i,]
    if(sum(yvec)==0)
    {loglik=loglik-lam[i]*(1-exp(-h*Tmax*J))}
    if(sum(yvec)>0)
    {# sort on yvec and t1
      t1vec=t1mat[i,]
      yt=cbind(yvec,t1vec)
      yt=yt[order(yt[,1],decreasing=FALSE),]
      yvec=yt[,1]
      t1vec=yt[,2]
      yij=yvec[-J]
      ymax=yvec[J]
      # terms for y[i,j] > 0 in site i
      for(j in 1:J)
      {
        y=yvec[j]
        if(y > 0)
        {term1=log(h)-h*t1vec[j]
        term2=(y-1)*log((exp(-h*t1vec[j])-exp(-h*Tmax)))
        term3=-sum(log(factorial(y-1)))
        loglik=loglik+term1+term2+term3
        }
      }
      #terms for all J equiv to sum over n
      term1=ymax*log(lam[i])
      term2=-lam[i]-h*Tmax*sum(ymax-yvec)
      term3=(J-1)*log(factorial(ymax))-sum(log(factorial(ymax-yij)))
      theta=lam[i]*exp(-h*Tmax*J)
      hgfun=genhypergeo(rep(ymax,(J-1)) + 1,ymax - yij + 1,theta)
      term4=log(hgfun)
      loglik=loglik+term1+term2+term3+term4}
  }
  return(-loglik)
}
