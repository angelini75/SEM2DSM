rm(list=ls()[])

# load lavaan model 
setwd("~/big/SEM2DSM1/Paper_4/data")
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
#z <- scale(z)  # make sure the column means are zero

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
alpha <- 0.5
# Range (a) ####
a <- 0.5 #e+05

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

# KronM function: adaptation of Kronecker function 
source("~/big/SEM2DSM1/Paper_4/kronM.R")

# next our approach with alpha != 0
MLIST$alpha <- alpha
MLIST$a <- a

z.all <- as.vector(t(z))  # compile to one big vector
RHO <- get.RHO(MLIST,h)
#SIGMA.all <- kronecker(RHO, SIGMA0)  # create covariance matrix of z.all
SIGMA.all <- kronM(RHO = RHO,RHO.I = diag(153),SIGMA0 = SIGMA0, sp = 1:9)
plotMat(SIGMA.all[(1+90):(90+90),(1+90):(90+90)])

L.all = chol(SIGMA.all)
logdetSIGMA.all = 2*sum(log(diag(L.all)))
SIGMA.all.inv <- chol2inv(L.all)

loglik.all <- -1/2*p*N*log(2*pi) - 1/2*logdetSIGMA.all - 
  1/2*t(z.all)%*%SIGMA.all.inv%*%z.all
loglik.all <- as.numeric(loglik.all)
# it runs without problems!

# now, let us define the objective function
objective_ML <- function(x, MLIST = MLIST) {
  MLIST <- x2MLIST(x = x, MLIST = MLIST)
  # compute Sigma.hat
  SIGMA0 <- computeSigmaHat.LISREL(MLIST = MLIST)
  RHO <- get.RHO(MLIST,h)
  SIGMA.all <- kronM(RHO = RHO,RHO.I = diag(153),SIGMA0 = SIGMA0, sp = 1:9)
    if (all(eigen(SIGMA.all)$values >0)) {
      L.all = chol(SIGMA.all)
    logdetSIGMA.all = 2*sum(log(diag(L.all)))
    SIGMA.all.inv <- chol2inv(L.all)
    objective <- -1 * (-1/2*p*N*log(2*pi) - 1/2*logdetSIGMA.all - 
      1/2*t(z.all)%*%SIGMA.all.inv%*%z.all)
    cat("objective = ", objective, "\n")
    objective
  } else {
    objective <- Inf
    objective
  }
}
# get the 51 starting values to feed x
lav.est <- parTable(my.fit.lv.ML)$est[parTable(my.fit.lv.ML)$free > 0]
start.x <- c(lav.est, alpha, a)

# optimizer of objective funtion 
sp.out  <- nlminb(start = start.x, objective = objective_ML, 
                  MLIST = MLIST, control = list(iter.max = 500, trace = 1))

round((start.x - sp.out$par),4)
x2MLIST(lav.out$par, MLIST)

# accuracy
MLIST.sp <- x2MLIST(sp.out$par, MLIST)

get.res <- function (m = NULL){
  A <- m$beta[1:9,10:18]
  B <- m$beta[1:9,1:9]
  I <- diag(nrow = 9, ncol = 9)
  IB.inv <- solve(I - B)
  sp <- s[,1:9]
  p <- s[,10:18]
  res <- matrix(data = NA, nrow = 153, ncol = 9)
  for(i in seq_along(p[,1])){
    res[i,] <- t(sp[i,] - (IB.inv %*% A %*% p[i,]))
  }
  colnames(res) <- colnames(sp)
  res
}
res <- get.res(m= MLIS.out)

rmse <- function(x){
  x = x^2
  y = sapply(as.data.frame(x), mean)
  z = sapply(y, sqrt)
  z
}

E <- data.frame(out1 = NA, lavaan = NA)

E[1:9,1] <- sapply(as.data.frame(res), FUN = rmse)
E[1:9,2] <- sapply(as.data.frame(res.lavaan), FUN = rmse)
rownames(E) <- colnames(s)[1:9]
E

# differences

A <- lavMLIST$beta[1:9,10:18] - MLIST.sp$beta[1:9,10:18]
colnames(A) <- colnames(s)[10:18]
rownames(A) <- colnames(s)[1:9]
levelplot(round(A,2), scale=list(x=list(rot=45)), 
          main = expression(Gamma))

B <- lavMLIST$beta[1:9,1:9] - MLIST.sp$beta[1:9,1:9]
colnames(B) <- colnames(s)[1:9]
rownames(B) <- colnames(s)[1:9]
levelplot(round(B,2), scale=list(x=list(rot=45)), 
          main = expression(Beta))

PSI <- lavMLIST$psi[1:9,1:9] - MLIST.sp$psi[1:9,1:9]
colnames(PSI) <- colnames(s)[1:9]
rownames(PSI) <- colnames(s)[1:9]
psiPlot <- levelplot(round(PSI,2), scale=list(x=list(rot=45)), 
          main = expression(Psi))
