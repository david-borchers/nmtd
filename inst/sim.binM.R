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
simest = matrix(rep(NA,nsim*2),nrow=nsim)
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

