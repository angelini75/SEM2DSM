set.seed(12345)
N <- 153  # number of observation locations
p <- 9   # number of variables
SIGMA <- SIGMA0  # to adjust with Bollen notation
plotMat(SIGMA)

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
