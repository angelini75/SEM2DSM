rm(list=ls()[])

# load lavaan model 
setwd("~/Documents/SEM2DSM1/Paper_4/data")
load("env.for.gerard.RData")
# lavaan model
fit <- my.fit.lv.ML
# estimated parameters with lavaan
MLIST <- lavTech(fit, "est")

# compute SIGMA0 (18x18)
SIGMA0 <- computeSigmaHat.LISREL(MLIST)
plotMat(SIGMA0)

# initialise matrix with standardised observations
# rows are locations, columns variables
z <- matrix(0, nrow=N, ncol=p)
for (i in 1:N) {
  z[i,] <- rnorm(p,0,1)  # assign random values
}
z <- scale(z)  # make sure the column means are zero

# first the Bollen method
L = chol(SIGMA)
logdetSIGMA = 2*sum(log(diag(L)))
SIGMA.inv <- chol2inv(L)
S <- 1/N*t(z)%*%z  # variance-covariance matrix, works better than var(z)

loglik <- -1/2*p*N*log(2*pi) - 1/2*N*logdetSIGMA - 
  1/2*N*sum(diag(S%*%SIGMA.inv))  # sum(diag()) gives the trace

# next our approach
z.all <- as.vector(t(z))  # compile to one big vector
SIGMA.all <- kronecker(RHO, SIGMA)  # create covariance matrix of z.all

L.all = chol(SIGMA.all)
logdetSIGMA.all = 2*sum(log(diag(L.all)))
SIGMA.all.inv <- chol2inv(L.all)

loglik.all <- -1/2*p*N*log(2*pi) - 1/2*logdetSIGMA.all - 
  1/2*t(z.all)%*%SIGMA.all.inv%*%z.all
loglik.all <- as.numeric(loglik.all)

# compare
loglik; loglik.all
