# Simulations without covariates
# ==============================
Rsites=100  #number of sites
Jsites=5    #multiple visits
Tsearch=3 #maximum time

#true parameters
lamt = 2
ht=0.4620981
paramt=c(lamt,ht)
init.paramt=c(log(lamt),log(ht))

# data for Binary:M
set.seed(123) # for reproducibility
bnm=as.matrix(generate.binM(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
# optimize
fit.binM=optim(init.paramt,nll.binM,R=Rsites,J=Jsites,Tmax=Tsearch,dat=bnm)
estpar.binM=fit.binM$par
exp(estpar.binM); paramt # compare estimates and true parameters


# Do some simulations to check:
nsim = 1000
npar = 3
simest = matrix(rep(NA,nsim*npar),nrow=nsim)
for(i in 1:nsim) {
  bns=as.matrix(generate.binM(paramt,R=Rsites,J=Jsites,Tmax=Tsearch))
  simfit.binM=optim(init.paramt,nll.binM,R=Rsites,J=Jsites,Tmax=Tsearch,dat=bns)
  simest[i,] = simfit.binM$par
}
parest.mean = c(mean(exp(simest[,1])), mean(exp(simest[,2])))
parest.mean; paramt

# some plots
par(mfrow=c(1,2))
# results for lamt
signif(mean(exp(simest[,1])),4)
signif(sd(exp(simest[,1])),4)
summary(exp(simest[,1]))
hist(exp(simest[,1]),25,main="lamt")
abline(v=mean(exp(simest[,1])))
abline(v=lamt,col="red")
boxplot(exp(simest[,1]),main="h")
abline(h=lamt,col="red")

# results for ht
signif(mean(exp(simest[,2])),4)
signif(sd(exp(simest[,2])),4)
summary(exp(simest[,2]))
hist(exp(simest[,2]),25,main="h")
abline(v=mean(exp(simest[,2])))
abline(v=ht,col="red")
boxplot(exp(simest[,2]),main="lamt")
abline(h=ht,col="red")


# Simulations with covariates
# ===========================
Rsites=100. #number of sites
Jsites=5    #number of visits
Tsearch=3 #maximum time

#true parameters
b0t =0.660878  #true intercept for log(lambda) 
b1t =0.255413 #true slope for log(lambda)
ht=0.4620981  #rate parameter

# simulation of the covariate x over R sites 
# hence the site-dependent abundance for the true parameters b0t and b1t
set.seed(123) # for reprodicibility
x=rnorm(R,0,1)
lambda=exp(b0t+b1t*x)
mean(lambda)

# true probabilites
1-exp(-ht*Tmax)      #P(individual detected)
1-exp(-mean(lambda)) #P(site occupancy)
paramt=c(b0t, b1t,ht)

# data for Binary:M
bnm=as.matrix(generate.binMcov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))

init.paramt=c(b0t, b1t, log(ht))
nll = nll.binMcov(param=init.paramt, R=Rsites, J=Jsites, Tmax=Tsearch,dat=bns, covar=x)
nll

# optimize
fit.binMcov=optim(init.paramt,nll.binMcov,R=Rsites,J=Jsites,Tmax=Tsearch,dat=bns,covar=x)
estpar.binMcov=fit.binMcov$par
# compare estimates and true parameters
c(estpar.binMcov[1:2],exp(estpar.binMcov[3]))
paramt 



# Do some simulations to check:
nsim = 1000
npar = 3
simest = matrix(rep(NA,nsim*npar),nrow=nsim)
for(i in 1:nsim) {
  bnm=as.matrix(generate.binMcov(paramt,R=Rsites,J=Jsites,Tmax=Tsearch,covar=x))
  simfit.binMcov=optim(init.paramt,nll.binMcov,R=Rsites,J=Jsites,Tmax=Tsearch,dat=bnm,covar=x)
  simest[i,] = simfit.binMcov$par
}
parest.mean = c(mean(simest[,1]), mean(simest[,2]), mean(exp(simest[,3])))
parest.mean; paramt

# some plots
# ----------
# results for b0t
par(mfrow=c(1,2))
# results for lamt
signif(mean(simest[,1]),4)
signif(sd(simest[,1]),4)
summary(simest[,1])
hist(simest[,1],25,main="b0t")
abline(v=mean(simest[,1]))
abline(v=b0t,col="red")
boxplot(simest[,1],main="h")
abline(h=b0t,col="red")

# results for b2t
par(mfrow=c(1,2))
# results for lamt
signif(mean(simest[,2]),4)
signif(sd(simest[,2]),4)
summary(simest[,2])
hist(simest[,2],25,main="b1t")
abline(v=mean(simest[,2]))
abline(v=b1t,col="red")
boxplot(simest[,2],main="h")
abline(h=b1t,col="red")


# results for ht
signif(mean(exp(simest[,3])),4)
signif(sd(exp(simest[,3])),4)
summary(exp(simest[,3]))
hist(exp(simest[,3]),25,main="h")
abline(v=mean(exp(simest[,3])))
abline(v=ht,col="red")
boxplot(exp(simest[,3]),main="lamt")
abline(h=ht,col="red")


# results for lambdas
nx = 100
xs = seq(min(x),max(x),length=nx)
lambda.est = matrix(rep(NA,nsim*nx),nrow=nsim)
for(i in 1:nsim) lambda.est[i,] = exp(simest[i,1] + simest[i,2]*xs)
mean.lambda.est = apply(lambda.est,2,mean)
lcl.lambda = apply(lambda.est,2,quantile,probs=0.025)
ucl.lambda = apply(lambda.est,2,quantile,probs=0.975)
ylim = range(lambda.est)
plot(xs,lambda.est[1,],type="l",col="gray",xlab="Covariate",ylab="Lambda",ylim=ylim)
for(i in 2:nsim) lines(xs,lambda.est[i,],col="gray")
lines(xs,mean.lambda.est)
lines(xs,lcl.lambda,lty=2)
lines(xs,ucl.lambda,lty=2)
lambda.true = exp(b0t + b1t*xs)
lines(xs,lambda.true,col="red")
