# generates binomial count data for model M3m with lambda constant
gM3Bmgen=function(param)
{lam=param[1]
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

# evaluates the log-likelihood for model M3m with lambda constant
gM3Bm=function(param)
{
  lam=exp(param[1])
  h=exp(param[2])
  loglik=0
  p=1-exp(-h*Tmax)
  ymat=gdat
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


# generates binomial count data for model M3m with lambda site-dependent
xM3Bmgen=function(param)
{b0=param[1]
b1=param[2]
lam=exp(b0+b1*vegHt)
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

# evaluates the log-likelihood for model M3m with lambda site-dependent
xgM3Bm=function(param)
{b0=param[1]
b1=param[2]
h=exp(param[3])
lam=exp(b0+b1*vegHt)
loglik=0
ymat=gdatx
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