Rsites=100 #number of sites
Tsearch=3 #maximum time

#true parameters
b0t =0.660878  #true intercept for log(lambda) 
b1t =0.255413 #true slope for log(lambda)
ht=0.4620981  #rate parameter

# simulation of the covariate vegHt over R sites 
# hence the site-dependent abundance for the true parameters b0t and b1t
set.seed(123) # for reprodicibility
x=rnorm(R,0,1)
lambda=exp(b0t+b1t*x)
mean(lambda)

# true probabilites
1-exp(-ht*Tmax)      #P(individual detected)
1-exp(-mean(lambda)) #P(site occupancy)
paramt=c(b0t, b1t,ht)

# data for Binary:S
bns=as.matrix(generate.binScov(paramt,R=Rsites,Tmax=Tsearch,covar=x))

init.paramt=c(b0t, b1t, log(ht))
nll = nll.binScov(param=init.paramt, R=Rsites, Tmax=Tsearch,dat=bns, covar=x)
nll
fit.binScov=optim(init.paramt,nll.binScov,R=Rsites,Tmax=Tsearch,dat=bns,covar=x)
estpar.binScov=fit.binScov$par
c(estpar.binScov[1:2],exp(estpar.binScov[3])); paramt # compare estimates and true parameters

# Do some simulations to check:
nsim = 1000
simest = matrix(rep(NA,nsim*3),nrow=nsim)
for(i in 1:nsim) {
  bns=as.matrix(generate.binScov(paramt,R=Rsites,Tmax=Tsearch,covar=x))
  simfit.binScov=optim(init.paramt,nll.binScov,R=Rsites,Tmax=Tsearch,dat=bns,covar=x)
  simest[i,] = simfit.binScov$par
}
parest.mean = c(mean(simest[,1]), mean(simest[,2]), mean(exp(simest[,3])))
parest.mean; paramt

# some plots
par(mfrow=c(1,2))
# results for b0
b0t
signif(mean(simest[,1]),4)
signif(sd(simest[,1]),4)
summary(simest[,1])
hist(simest[,1],25,main="b0")
abline(v=mean(simest[,1]))
abline(v=b0t,col="red")
boxplot(simest[,1],main="b0")
abline(h=b0t,col="red")

# results for b1
b1t
signif(mean(simest[,2]),4)
signif(sd(simest[,2]),4)
summary(simest[,2])
hist(simest[,2],25,main="b1")
abline(v=mean(simest[,2]))
abline(v=b1t,col="red")
boxplot(simest[,2],main="b1")
abline(h=b1t,col="red")

# results for ht
ht
signif(mean(exp(simest[,3])),4)
signif(sd(exp(simest[,3])),4)
summary(exp(simest[,3]))
hist(exp(simest[,3]),25,main="h")
abline(v=mean(exp(simest[,3])))
abline(v=ht,col="red")
boxplot(exp(simest[,3]),main="h")
abline(h=ht,col="red")

