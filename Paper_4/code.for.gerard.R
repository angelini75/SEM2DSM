# library(lavaan)
# HS.model <- ' visual  =~ x1 + x2 + x3
# textual =~ x4 + x5 + x6
# speed   =~ x7 + x8 + x9 
# textual ~ speed + visual
# speed ~ visual'
# 
# fit <- sem(HS.model, data=HolzingerSwineford1939, estimator = "ML") # or ML
library(lavaan)
library(lattice)
library(sp)
library(gstat)
rm(list=ls()[])

# load lavaan model 
setwd("~/big/SEM2DSM1/Paper_4/data")
load("env.for.gerard.RData")
fit <- my.fit.lv.ML
# initial list of model matrices
MLIST <- lavTech(fit, "start")

# sample covariance matrix from lavaan. It uses N
S     <- fit@SampleStats@cov[[1]]
S[1:3,1:3]
# similar to var-covar matrix of observation when using N-1
cov(s)[1:3,1:3] # s are observed standardized data
S[1:3,1:3] - cov(s)[1:3,1:3]

# RHO is an identity matrix of N x N (for alpha = 0) 
RHO <- diag(153)
# Lets create  a sample var-covar matrix of (q+p)N x (q+p)N
Si <- kronecker(RHO,S)
plotMat(Si[1:72,1:72])

# inverse and logdet of Si. 
cS <- chol(Si)
S.inv <- chol2inv(cS)
diag.cS <- diag(cS)
S.logdet <- sum(log(diag.cS * diag.cS))
# objective function 'ML'

# give estimates as starting values
start.x <- parTable(fit)$est[ parTable(fit)$free > 0 ]
x <- start.x

# MLIST has all the (free and fixed) parameters needed to estimate SIGMA
# and x2SIGMA is a function to insert the starting values in MLIST
MLIST <- x2MLIST(x = x, MLIST = MLIST)

# compute SIGMA0 and SIGMA
SIGMA0 <- computeSigmaHat.LISREL(MLIST = MLIST) # 18x18
plotMat(SIGMA0)
SIGMA <- kronecker(RHO, SIGMA0)
plotMat(SIGMA[1:72,1:72])

# function used by lavaan + RHO
objective_ML <- function(x, MLIST) {
  MLIST <- x2MLIST(x = x, MLIST = MLIST)
  # compute Sigma.hat
  SIGMA0 <- computeSigmaHat.LISREL(MLIST = MLIST)
  if (all(eigen(SIGMA0)$values >0)) {
    SIGMA <- kronecker(RHO, SIGMA0)
    nvar <- NROW(SIGMA)
    cS <- chol(SIGMA)
    Sigma.inv <- chol2inv(cS)
    diag.cS <- diag(cS)
    logdet <- sum(log(diag.cS * diag.cS))
    objective <- logdet + sum(diag(Si %*% Sigma.inv)) - S.logdet - nvar
    cat("objective = ", objective, "\n")
    objective
  } else {
    objective <- Inf
    objective
  }
}

# estimate parameters (it takes about two hours)
# You need to run 150 iteration (eval.max = 150) to get small differences
# between lavaan and lavaan+RHO
out.ML  <- nlminb(start = start.x, objective = objective_ML, 
                  MLIST = MLIST, control = list(eval.max = 150, trace = 1))
# Differences between lavaan and lavaan+RHO
round((parTable(my.fit.lv.ML)$est[parTable(my.fit.lv.ML)$free > 0] - out.ML$par[1:51]),4)
MLIST.out <- x2MLIST(out.ML$par, MLIST)

