rm(list=ls()[])
library(lavaan)

# load lavaan model 
setwd("~/big/SEM2DSM1/Paper_4/data")
load("env.for.gerard.RData")
load("lavaan.model.RData")
ks <- read.csv("ks.csv")
ks <- ks[,colnames(s)]

# lavaan model
fit <- my.fit.lv.ML

# estimated parameters with lavaan
MLIST <- lavTech(fit, "est")

# compute SIGMA0 (18x18)
SIGMA0 <- computeSigmaHat.LISREL(MLIST)
plotMat(SIGMA0)

# initialise matrix with standardised observations
# rows are locations, columns variables
z <- as.matrix(ks)

################################################################################
ks <- read.csv("ks.csv")[,-1] # standardized data
ST <- read.csv("STt.ks.csv")
ks$Y2 <- ks$Y * ST$std.dev[21] / ST$std.dev[20]

coordinates(ks) <- ~X+Y2
# h = n x n matrix of distances between samples 
h <- spDists(ks)

# Define the parameters alpha and range, and estimate distance matrix:
# Sill to nugget ratio (alpha) ####
# alpha = C/(C0 + C) = 0.8/(0.2 + 0.8) 
alpha <- 0.5
# Range (a) ####
a <- 0.7 #e+05

# Function to include free parameters into the matrix structure MLIST
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
    RHO[i] <- (1-alpha) * exp(-h[i]/a)
  }
  diag(RHO) <- 1
  RHO
}

# KronM function: adaptation of Kronecker function [NOT NECESSARY]
#source("~/big/SEM2DSM1/Paper_4/kronM.R")

# next our approach with alpha != 0
MLIST$alpha <- alpha
MLIST$a <- a

z.all <- as.vector(z)  # compile to one big vector
RHO <- get.RHO(MLIST,h)
SIGMA.all <- kronecker(SIGMA0, RHO)  # create covariance matrix of z.all
plotMat(SIGMA.all[1000:1700,1000:1700])
# SIGMA.all2 <- kronM(RHO = RHO,RHO.I = RHO,SIGMA0 = SIGMA0, sp = 1:9, cov = 10:18)
# plotMat(SIGMA.all1[1000:1500,1000:1500]-SIGMA.all2[1000:1500,1000:1500])


# Choleschy decomposition
L.all = chol(SIGMA.all)
logdetSIGMA.all = 2*sum(log(diag(L.all)))
SIGMA.all.inv <- chol2inv(L.all)

p=18
N=153

# to define the objective function:
objective_ML <- function(x, MLIST = MLIST) {
  MLIST <- x2MLIST(x = x, MLIST = MLIST)
  # compute Sigma.hat
  SIGMA0 <- computeSigmaHat.LISREL(MLIST = MLIST)
  RHO <- get.RHO(MLIST,h)
  if (all(eigen(SIGMA0)$values >0) & (all(eigen(RHO)$values >0))) {
    SIGMA.all <- kronecker(SIGMA0, RHO)
    #SIGMA.all <- kronM(RHO = RHO,RHO.I = RHO,SIGMA0 = SIGMA0, sp = 1:9, cov = 10:18)
    #  if (all(eigen(SIGMA.all)$values >0)) {
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

start.x <- c(lav.est, 0.8, 0.5)
sp.ou <- nlminb(start = start.x, objective = objective_ML, 
                MLIST = MLIST, control = list(iter.max = 500, trace = 1))

round((start.x - sp.ou$par),4)
MLIST.out <- x2MLIST(sp.ou$par, MLIST)






