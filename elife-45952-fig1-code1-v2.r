# Figure 1: DFEs
##### DFE (all) #####

data <- gtsB

hist(data$mean.w[data$effect!="S"], breaks=seq(0.6,1.1, by=0.02), 
     xlab="Relative fitness", las=1, col=rgb(0,0,1,1/4),
     main="", ylab="", ylim = c(0,14))
hist(data$mean.w[data$effect=="S"], breaks=seq(0.6,1.1, by=0.02), col=rgb(1,0,0,1/4), add=TRUE)
rug(data$mean.w[data$effect=="STOP"], col="grey40", lwd=2)
abline(v=1.004, lty = 2)  # mean fitness of SBW25 competed against itself from all assays
abline(v=0.9919863561, lty=3) #lower 95% CI SBW25
abline(v=1.017478972, lty=3) #upper 95% CI SBW25


## re-sample and calculate a KS statistic (KS assumes continuous ditribution)
NS.counts = hist(NS.data$mean.w, breaks=seq(from=0.65, to=1.15, by=0.05), plot=F)$counts
mids = hist(NS.data$mean.w, breaks=seq(from=0.65, to=1.15, by=0.05), plot=F)$mids

N.permut = 10000
####### below: bootstrapping the data
KS.D.dist = vector(length = N.permut)
for (n in 1:N.permut) {
  NS.sample = sample(x = NS.data$mean.w, size = Num.S, replace = T)
  ### comparing 2 distns
  ### and S is smaller so to make sample sizes equal take S size
  # hist(NS.sample, breaks=seq(from=0.5, to=1.2, by=0.05), plot=T)
  NS.sample.counts = hist(NS.sample, breaks=seq(from=0.65, to=1.15, by=0.05), plot=F)$counts
  M = as.table(rbind(NS.counts,NS.sample.counts))
  
  KS.D.dist[n] = ks.test(x = NS.data$mean.w, y = NS.sample)$statistic
}
KS.D.real = ks.test(x = NS.data$mean.w, y = S.data$mean.w)$statistic
(KS.P.value = sum(KS.D.real<KS.D.dist)/N.permut)


## K-S test for beneficial mutations
N.permut = 10000
KS.D.dist = vector(length = N.permut)
for (n in 1:N.permut) {
  NS.sample = sample(x = NS.data$mean.w[NS.data$w.effect==1], size = length(S.data$mean.w[S.data$w.effect==1]), replace = T)
  KS.D.dist[n] = ks.test(x = NS.data$mean.w[NS.data$w.effect==1], y = NS.sample)$statistic
}
KS.D.real = ks.test(x = NS.data$mean.w[NS.data$w.effect==1], y = S.data$mean.w[S.data$w.effect==1])$statistic
(KS.P.value = sum(KS.D.real<KS.D.dist)/N.permut)

## Compare N and NS beneficial distributions ####
## plot smoothed distributions
plot(N.dens<-density(NS.data$mean.w[NS.data$w.effect==1], bw=0.012, 
                     from = 1, to = 1.15), ylim=c(0,22), 
     xlim=c(1,1.15), lwd=2, main = "DFE of beneficial mutations", xlab = "Relative fitness")
lines(density(S.data$mean.w[S.data$w.effect==1], bw = N.dens$bw, from = 1, to = 1.15), col="red", lwd=2)
legend(1.1, 20, legend = c("Non-synonymous", "Synonymous"), col = c(1,2), lwd = 2, bty = "n")
text(1.11, 5, 
     "Permutation test of K-S statistics: \n distributions are not significantly different \n(p = 0.634)",
     pos = 3)

## Testing the shape of the distribution of beneficial mutation ####
## Use both the synonymous and non-synonymous together because their beneficial distributions are not significantly different.
## Define the mutations of interest:
s = s.ben # syn and non-synonymous

## Instead of measuring fitness relative to the wild type, 
## we shift the distribution relative to the smallest observed 
## selection coefficient (method from Beisel et al 2007)
X = s-min(s)
X = X[X!=0]
X = sort(X)

LL.GPD <- function(par){ # equations from Beisel et al 2007
  # X is the adjusted selection coefficients, i.e. the data
  X = X
  # tau is a parameter of the GPD
  tau = par[1]
  # kappa is a parameter of the GDP
  kappa = par[2]
  # n_1 is the number of observations
  n_1 = length(X)
  LL.a = -(n_1)*log(tau)
  LL.b = NULL
  # for kappa>0
  if(kappa>0) {
    LL.b = -(kappa+1)/kappa * sum (log(1+(kappa*X)/tau))
  }
  # for kappa<0
  if(kappa<0) {
    if(sum(X>-tau/kappa)){LL.b = -1000000}
    else{LL.b = -(kappa+1)/kappa * sum (log(1+(kappa*X)/tau))}
  }
  # for kappa=0
  if(kappa==0) {
    LL.b = -(1/tau)*sum(X)
  }
  LL = LL.a + LL.b
  return(-LL)
}

## Special case kappa=0 (exponential distribution)
LL.Exp <- function(par){ # equations from Beisel et al 2007
  # X is the adjusted selection coefficients, i.e. the data
  X = X
  # tau is a parameter of the GPD
  tau = par
  # kappa is a parameter of the GDP
  kappa = 0
  
  # n_1 is the number of observations
  n_1 = length(X)
  LL.a = -(n_1)*log(tau)
  LL.b = -(1/tau)*sum(X)
  LL = LL.a + LL.b
  return(-LL)
}


#1# Optimization ####
# install.packages("GenSA")
require(GenSA)

start.par = c(.1, 0) # starting parameter values for the optimization
## 1) Find tau, When kappa = 0
Exp.opt = GenSA(par = start.par[1], fn = LL.Exp, lower = 0.000001, upper = 100)
(Exp.opt.tau = Exp.opt$par)
(LL.Exp.opt = -Exp.opt$value)

## 2) Find best tau and kappa
GPD.opt = GenSA(par = start.par, fn = LL.GPD, lower = c(0.000001, -100), upper = c(100, 100))
(GPD.opt.tau = GPD.opt$par[1])
(GPD.opt.kappa = GPD.opt$par[2])
(LL.GPD.opt = -GPD.opt$value)

## Likelihood ratio ####
(LRT = LL.Exp.opt/LL.GPD.opt)

## Find a null distribution of likelihood ratio ####
B = 10000
LR.null.dist = vector(length = B)
for(b in 1:B){
  X = rexp(n = length(X), rate = 1/Exp.opt.tau)
  Exp.null.opt = optim(par = start.par[1], fn = LL.Exp, method = "Brent", lower = 0.000001, upper = 100)
  LL.Exp.null.opt = -Exp.null.opt$value
  ## 2) Find best tau and kappa
  GPD.null.opt = optim(par = c(Exp.opt.tau, 0), fn = LL.GPD, method = "L-BFGS-B", lower = c(0.000001, -100), upper = c(100, 100))
  # print(b)
  # GPD.null.opt = GenSA(par = start.par, fn = LL.GPD, lower = c(0.000001, -100), upper = c(100, 100))
  LL.GPD.null.opt = -GPD.null.opt$value
  ## Likelihood ratio
  LR.null.dist[b] = LL.Exp.null.opt/LL.GPD.null.opt
}
(P.value = sum(LRT>LR.null.dist)/B) 


## Plot data and best fit distribution ####
install.packages("fExtremes")
require(fExtremes)

## Results plot ####
hist(X, freq = F, main = 'All beneficial mutations', breaks = 15, ylim = c(0, 35), 
     xlab = "s", col = c("mediumpurple1"))
y = seq(from=0.001, to=0.15, by=0.001)
lines(y, dgpd(x = y, xi = GPD.opt.kappa, mu = 0, beta = GPD.opt.tau), col='black', lwd=2)
#lines(y, dexp(x = y, rate = 1/Exp.opt.tau), col='black', lwd=2)
legend(0.04, 25, legend = c( 
  expression(paste("Generalized Pareto (", kappa, " is fit)"))), col = c("black"), lwd = 2, bty = "n")
text(0.08, 15, labels = expression(paste("Weibull (", kappa, " = -0.370, ", tau, " = 0.0330)")))

