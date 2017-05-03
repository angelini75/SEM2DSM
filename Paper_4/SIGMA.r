library(lavaan)
library(lattice)
library(sp)
library(gstat)


rm(list=ls()[ls()!="my.fit.lv.ML"])
# utils
name <- function(x) { as.data.frame(names(x))}
plotMat <- function(matrix=matrix){
  plotMat <- t(apply(matrix,2,rev)) # rotate 90 degrees clockwise
  library(lattice)
  print(levelplot(plotMat, ylab="row", xlab="column"))
}
################################################################################
################ HOW LAVAAN ESTIMATES THE COEFICIENTS? #########################
################################################################################
# load lavaan model 
setwd("~/Documents/SEM2DSM1/Paper_4/data")
load("lavaan.model.RData")

# adapted from https://groups.google.com/forum/#!searchin/lavaan/lavaan$20function$20github%7Csort:relevance/lavaan/0Hitqi3k_do/pYfyoocABwAJ

# This are the estimates of our model
lavTech(my.fit.lv.ML, "est")

# lavMLIST "start" have the starting values of free parameters
lavMLIST <- lavTech(my.fit.lv.ML, "start")

# This is the sample covariance matrix called S
S <- my.fit.lv.ML@SampleStats@cov[[1]]

# inverse and log(det(S))
cS <- chol(S)
S.inv <- chol2inv(cS)
diag.cS <- diag(cS)
S.logdet <- sum(log(diag.cS * diag.cS))

# x2MLIST (renamed as x2lavMLIST) fill in new values (from x) in the 
# model matrices (lavMLIST)
par.list <- inspect(my.fit.lv.ML) # the same than lavMLIST
foo <- c(101:151) #51 parameters
x2lavMLIST <- function(x, lavMLIST = lavMLIST) {
  lambda.x <- x[as.vector(par.list$lambda)[as.vector(par.list$lambda)!=0]]
  theta.x <- x[as.vector(par.list$theta)[as.vector(par.list$theta)!=0]]
  psi.x <- x[as.vector(par.list$psi)[as.vector(par.list$psi)!=0]]
  beta.x <- x[as.vector(par.list$beta)[as.vector(par.list$beta)!=0]]
  
  lavMLIST$lambda[which(as.vector(par.list$lambda)!=0)]     <- lambda.x
  lavMLIST$theta[which(as.vector(par.list$theta)!=0)]       <- theta.x
  lavMLIST$psi[which(as.vector(par.list$psi)!=0)]           <- psi.x
  lavMLIST$beta[which(as.vector(par.list$beta)!=0)]           <- beta.x
  lavMLIST
}
x2lavMLIST(foo, lavMLIST = lavMLIST) 

# function to get inverse(I - B) = solve(I - B). We leave it as it is.
# (I think it solves issues for non-positive definite -PD- matrices)
get.IB.inv <- function (lavMLIST = NULL) {
  BETA <- lavMLIST$beta
  nr <- nrow(lavMLIST$psi)
  if (!is.null(BETA)) {
    tmp <- -BETA
    tmp[lav_matrix_diag_idx(nr)] <- 1
    IB.inv <- solve(tmp)
  }
  else {
    IB.inv <- diag(nr)
  }
  IB.inv
}

# Since lavaan uses Sigma-hat intead of residuals in the ML function,
# it needs a function to estimate it
computeSigmaHat.LISREL <- function (lavMLIST = NULL, delta = TRUE) 
{
  LAMBDA <- lavMLIST$lambda # links between observed var and latent variables
  nvar <- nrow(LAMBDA) 
  PSI <- lavMLIST$psi # this is PSI and PHI together
  THETA <- lavMLIST$theta # measurement errors of both X's and Y's
  BETA <- lavMLIST$beta # GAMMA and BETA coeficients
  library(lavaan)
  if (is.null(BETA)) {
    LAMBDA..IB.inv <- LAMBDA
  }
  else { 
    lavaan:::.internal_get_IB.inv(MLIST = lavMLIST)
    IB.inv <- get.IB.inv(lavMLIST = lavMLIST)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv 
  }
  VYx <- tcrossprod(LAMBDA..IB.inv %*% PSI, LAMBDA..IB.inv) + THETA # key function
  if (delta && !is.null(lavMLIST$delta)) {
    DELTA <- diag(lavMLIST$delta[, 1L], nrow = nvar, ncol = nvar)
    VYx <- DELTA %*% VYx %*% DELTA
  }
  VYx
}

# objective function 'ML'

objective_ML <- function(x, lavMLIST = lavMLIST) {
  lavMLIST <- x2lavMLIST(x = x, lavMLIST = lavMLIST)
  # compute Sigma.hat
  Sigma <- computeSigmaHat.LISREL(lavMLIST = lavMLIST)
  if (all(eigen(Sigma)$values >0)) {
    nvar <- NROW(Sigma)
    cS <- chol(Sigma)
    Sigma.inv <- chol2inv(cS)
    diag.cS <- diag(cS)
    logdet <- sum(log(diag.cS * diag.cS))
    objective <- logdet + sum(diag(S %*% Sigma.inv)) - S.logdet - nvar
    cat("objective = ", objective, "\n")
    objective
  } else {
    objective <- Inf
    objective
  }
}
# get the 51 starting values to feed x
(lav.start.x <- parTable(my.fit.lv.ML)$start[parTable(my.fit.lv.ML)$free > 0])

# optimizer of objective funtion 
lav.out  <- nlminb(start = lav.start.x, objective = objective_ML, 
                  lavMLIST = lavMLIST)
# Now, we assign the coefficients to the lavMLIST matrices
lavMLIST <- x2lavMLIST(lav.out$par, lavMLIST)


################################################################################
########## HOW TO ESTIMATES THE COEFICIENTS WITH SPATIAL PROCESS? ##############
################################################################################

# Load data with locations in the same order that is defined in lavaan
sp <- rownames(inspect(my.fit.lv.ML, "est")$lambda)[1:9] #names of soil properties
covar <- rownames(inspect(my.fit.lv.ML, "est")$lambda)[10:18] # names of covariates
ks <- read.csv("ks.csv")[,-1] # standardized data
s <- as.matrix(ks[,sp]) # soil properties
p <- as.matrix(ks[,covar]) #covariates


# Define the parameters alpha and range, and estimate distance matrix:
# Sill to nugget ratio (alpha) ####
# alpha = C/(C0 + C) = 0.8/(0.2 + 0.8) 
alpha <- 0 
# Range (a) ####
a <- 1 #e+05

# Distance (h) ####
# make ks spatial objetc
coordinates(ks) <- ~X+Y
# h = n x n matrix of distances between samples 
h <- spDists(ks) # in reality this will vary 

# GAMMA, BETA, PSI and THETA matrices ####
# We use lavaan estimates ("est") as starting values 
# Matrix of BETA coefficients
B <- inspect(my.fit.lv.ML, "est")$beta[1:9,1:9] # s x s
# Identity matrix
I <- diag(nrow = 9, ncol = 9)
# Matrix of GAMMA coefficients
A <- inspect(my.fit.lv.ML, "est")$beta[1:9,10:18] # s x p
# Matrix of Psi coefficients (model error variance-covariance)
PSI <- inspect(my.fit.lv.ML, "est")$psi[1:9,1:9] # s x s
# Matrix of measurement error (Epsylon)
TH <- inspect(my.fit.lv.ML, "est")$theta[1:9,1:9] # s x s
# (I - B)^{inverse}
get.IB.inv <- function (m = NULL) {
  BETA <- m$B
  nr <- nrow(m$PSI)
  if (!is.null(BETA)) {
    tmp <- -BETA
    tmp[lav_matrix_diag_idx(nr)] <- 1
    IB.inv <- solve(tmp)
  }
  else {
    IB.inv <- diag(nr)
  }
  IB.inv
}

tmp <- -B
tmp[lavaan:::lav_matrix_diag_idx(nrow(PSI))] <- 1
IB.inv <- solve(tmp)

# Now, we create matrices with free and fixed parameters
# MLIST ####
MLIST <- list(A = A, B = B, PSI = PSI, alpha = alpha, a = a, TH = TH,
              sp = s, covar = p, h = h)


#### FUNTIONS ####
# To estimate the parameters we need four matrices:
# SIGMA0, that is the s x s var-covar matrix at h=0
# RHO, based on the variogram model
# Residuals (res), an sn vector of observed - predicted
# SIGMA, that is the sn x sn var-covar matrix at h=h

# Get SIGMA0
# 
get.SIGMA0 <- function (m = NULL) {
  PSI <- m$PSI
  THETA <- m$TH
  IB.inv <- get.IB.inv(m = MLIST)
  IB.PSI.IBinv <- tcrossprod(IB.inv %*% PSI, IB.inv) + THETA
  IB.PSI.IBinv
}
plotMat(get.SIGMA0(m = MLIST))

# Get RHO 
get.RHO <- function(m = NULL) {
  a <- m$a
  n <- nrow(m$h)
  h <- m$h
  alpha <- m$alpha
  RHO <- matrix(rep(NA,n^2), nrow = n)
  for(i in seq_along(RHO)) {
    RHO[i] <- alpha * exp(-h[i]/a)
  }
  diag(RHO) <- 1
  RHO
}
plotMat(get.RHO(m = MLIST)[1:27,1:27])

# Get residuals 
get.res <- function (m = NULL){
  A <- m$A
  B <- m$B
  I <- diag(nrow = 9, ncol = 9)
  IB.inv <- get.IB.inv(m = m)
  s <- m$sp
  p <- m$covar
  res <- matrix(data = NA, nrow = 153, ncol = 9)
  for(i in seq_along(p[,1])){
    res[i,] <- t(s[i,] - (IB.inv %*% A %*% p[i,]))
  }
  colnames(res) <- colnames(m$sp)
  res
}
head(get.res(m = MLIST)) # Warning: it returns a data frame

# Get SIGMA 
# SIGMA is just the kronecker product of RHO and SIGMA0
plotMat(kronecker(get.RHO(MLIST),
                  get.SIGMA0(MLIST))[1:60,1:60])

# Now, a function to fill in MLIST with new values
new.par.list <- list(A = par.list$beta[1:9,10:18], 
                     B = par.list$beta[1:9,1:9],
                     PSI = par.list$psi[1:9,1:9],
                     alpha = 52,
                     a = 53)
x2MLIST <- function(x, MLIST) {
  par <- list(A = x[1:25],
            B = x[26:37],
            PSI = x[38:51])#,
#            alpha = x[52],
#            a = x[53])
  MLIST$A[which(as.vector(new.par.list$A)!=0)]          <- par$A
  MLIST$B[which(as.vector(new.par.list$B)!=0)]          <- par$B
  MLIST$PSI[which(as.vector(new.par.list$PSI)!=0)]      <- 
    par$PSI[c(1,2,3,4, 2 ,5,6, 3 , 6 ,7,8,9,10,11,12, 4 , 11 ,13,14)]
  # MLIST$alpha[which(as.vector(new.par.list$alpha)!=0)]  <- par$alpha
  # MLIST$a[which(as.vector(new.par.list$a)!=0)]          <- par$a
  MLIST
}
# check how it works
MLIST2 <- x2MLIST(x = 1:51, MLIST)
MLIST2[1:3]
# MLIST2$A - inspect(my.fit.lv.ML, "est")$beta[1:9,10:18]
# MLIST2$B - inspect(my.fit.lv.ML, "est")$beta[1:9,1:9]
# MLIST2$PSI - inspect(my.fit.lv.ML, "est")$psi[1:9,1:9]

# objective function 'ML'
objective_ML <- function(x, MLIST) {
  MLIST <- x2MLIST(x = x, MLIST = MLIST)
  # compute Sigma.hat
  res <- get.res(MLIST = MLIST)
  SIGMA0 <- get.IB.PSI.IBinv(MLIST = MLIST)
  RHO <- get.RHO(MLIST = MLIST)
  SIGMA <- kronecker(RHO, SIGMA0)
  if (all(eigen(SIGMA0, only.values = T)$values > 0)) {
    nvar <- NROW(SIGMA)
    cS <- chol(SIGMA)
    Sigma.inv <- chol2inv(cS)
    diag.cS <- diag(cS)
    logdet <- sum(log(diag.cS * diag.cS))
    objective <- logdet + t(res) %*% Sigma.inv %*% res
    cat("objective = ", objective, "\n")
    objective
  } else {
    objective <- Inf
    objective
  }
}
# get starting values
start.x <- inspect(my.fit.lv.ML, "est")$beta[1:9,10:18][
  which(as.vector(new.par.list$A)!=0)]
start.x <- append(start.x, inspect(my.fit.lv.ML, "est")$beta[1:9,1:9][
  which(as.vector(new.par.list$B)!=0)])
starting.psi <- partable(my.fit.lv.ML)$est[
  unique(as.vector(lavTech(my.fit.lv.ML,"partable")$psi[1:9,1:9]))
]
start.x <- append(start.x, starting.psi)
#start.x <- append(start.x, MLIST$alpha[which(as.vector(new.par.list$alpha)!=0)])
#start.x <- append(start.x, MLIST$a[which(as.vector(new.par.list$a)!=0)] )

# estimate parameters ML
out.ML  <- nlminb(start = start.x, objective = objective_ML, 
                  MLIST = MLIST, control = list(eval.max = 5000, trace = 1))

x2MLIST(out.ML$par, MLIST)


det(SIGMA)

