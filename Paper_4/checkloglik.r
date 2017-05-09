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
z <- s
z <- scale(z)  # make sure the column means are zero

# first the Bollen method
N = 153
p = 18
L = chol(SIGMA0)
logdetSIGMA = 2*sum(log(diag(L)))
SIGMA.inv <- chol2inv(L)
S <- 1/N*t(z)%*%z  # variance-covariance matrix, works better than var(z)
lavS <- my.fit.lv.ML@SampleStats@cov[[1]]

S - lavS

loglik <- -1/2*p*N*log(2*pi) - 1/2*N*logdetSIGMA - 
  1/2*N*sum(diag(S%*%SIGMA.inv))  # sum(diag()) gives the trace

# next our approach
z.all <- as.vector(t(z))  # compile to one big vector
RHO <- diag(153)
SIGMA.all <- kronecker(RHO, SIGMA0)  # create covariance matrix of z.all

L.all = chol(SIGMA.all)
logdetSIGMA.all = 2*sum(log(diag(L.all)))
SIGMA.all.inv <- chol2inv(L.all)

loglik.all <- -1/2*p*N*log(2*pi) - 1/2*logdetSIGMA.all - 
  1/2*t(z.all)%*%SIGMA.all.inv%*%z.all
loglik.all <- as.numeric(loglik.all)

# compare
loglik; loglik.all

################################################################################
ks <- read.csv("ks.csv")[,-1] # standardized data
coordinates(ks) <- ~X+Y
# h = n x n matrix of distances between samples 
h <- spDists(ks)

# Define the parameters alpha and range, and estimate distance matrix:
# Sill to nugget ratio (alpha) ####
# alpha = C/(C0 + C) = 0.8/(0.2 + 0.8) 
alpha <- 1
# Range (a) ####
a <- 1 #e+05

# Change x2MLIST to add alpha and a parameters
x2MLIST <- function(x, MLIST) {
  lambda.x <- x[as.vector(par.list$lambda)[as.vector(par.list$lambda)!=0]]
  theta.x <- x[as.vector(par.list$theta)[as.vector(par.list$theta)!=0]]
  psi.x <- x[as.vector(par.list$psi)[as.vector(par.list$psi)!=0]]
  beta.x <- x[as.vector(par.list$beta)[as.vector(par.list$beta)!=0]]
  alpha.x <- x[52]
  a.x <- x[53]
  
  MLIST$lambda[which(as.vector(par.list$lambda)!=0)]     <- lambda.x
  MLIST$theta[which(as.vector(par.list$theta)!=0)]       <- theta.x
  MLIST$psi[which(as.vector(par.list$psi)!=0)]           <- psi.x
  MLIST$beta[which(as.vector(par.list$beta)!=0)]         <- beta.x
  MLIST$alpha                                            <- alpha.x
  MLIST$a                                                <- a.x
  MLIST
}

# get.RHO function
get.RHO <- function(MLIST = NULL, h = h) {
  a <- MLIST$a
  n <- nrow(h)
  alpha <- MLIST$alpha
  RHO <- matrix(rep(NA,n^2), nrow = n)
  for(i in seq_along(RHO)) {
    RHO[i] <- alpha * exp(-h[i]/a)
  }
  diag(RHO) <- 1
  RHO
}

# next our approach with alpha != 0
MLIST$alpha <- alpha
MLIST$a <- a

z.all <- as.vector(t(z))  # compile to one big vector
RHO <- get.RHO(MLIST,h)
SIGMA.all <- kronecker(RHO, SIGMA0)  # create covariance matrix of z.all
plotMat(SIGMA.all[1:72,1:72])

L.all = chol(SIGMA.all)
logdetSIGMA.all = 2*sum(log(diag(L.all)))
SIGMA.all.inv <- chol2inv(L.all)

loglik.all <- -1/2*p*N*log(2*pi) - 1/2*logdetSIGMA.all - 
  1/2*t(z.all)%*%SIGMA.all.inv%*%z.all
loglik.all <- as.numeric(loglik.all)
