library(lattice)

# 9 soil properties
m <- 9
# 9 covariates
q <- 9
# n = samples
n=153 

# Set alpha ####
# alpha = C/(C0 + C) = 0.8/(0.2 + 0.8) 
alpha <- 0.8 
# Set range ####
a <- 1e+05

# Get distances (h) from locations ####
# Load locations 
setwd("~/Documents/SEM2DSM1/Paper_4/data") # change folder if necessary
ks <- read.csv("ks.csv")[,-1] # standardized data
STt.ks <- read.csv("STt.ks.csv") # mean and sd of each variable
# function to unstandardize 
unstd<- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] * st[i,3]) + st[i,2]
  }
  y
}
# unstandardize coordinates X and Y
R <- unstd(x = ks[,c(20,21)],st = STt.ks[c(20,21),])

library(sp)
library(gstat)
coordinates(R) <- ~X+Y
#define crs
NAD83.KS.N <- CRS("+init=epsg:2796")
proj4string(R) <- NAD83.KS.N

# h = n x n matrix of distances between samples 
h <- spDists(R) # in reality this will vary 

# Get RHO ####
RHO <- matrix(rep(NA,n^2), nrow = n)
for(i in seq_along(RHO)) {
  RHO[i] <- alpha * exp(-h[i]/a)
}
diag(RHO) <- 1
RHO

# Get SIGMA0 ####
#from model 4
# initial list of model matrices
    # MLIST <- lavTech(my.fit.lv.ML, "start")
    # # sample covariance matrix
    # #S     <- fit@SampleStats@cov[[1]]
    # 
    # # inverse and logdet of S
    # cS <- chol(S)
    # S.inv <- chol2inv(cS)
    # diag.cS <- diag(cS)
    # S.logdet <- sum(log(diag.cS * diag.cS))

# x is the parameter vector: fill in new values in the model matrices
par.list <- inspect(my.fit.lv.ML)

# The starting values of Sigma zero (SIGMA0) 
# Matrix of Beta coefficients
B <- inspect(my.fit.lv.ML, "start")$beta[1:9,1:9]
# Identity matrix
I <- diag(nrow = 9, ncol = 9)
# Matrix of Gamma coefficients
#A <- lavTech(my.fit.lv.ML, "start")$beta[1:9,10:18]
A <- inspect(my.fit.lv.ML, "start")$beta[1:9,10:18] # 9 soil properties x 9 covariates
# Matrix of Psi coefficients (model error variance-covariance)
PSI <- inspect(my.fit.lv.ML, "start")$psi[1:9,1:9]
# Matrix of measurement error (Epsylon)
TH <- inspect(my.fit.lv.ML, "start")$theta[1:9,1:9]

IB <- solve(I - B)

# Get residuals ####
p <- as.matrix(ks[,colnames(A)]) #covariates
s <- as.matrix(ks[,names(ks[2:10])])
res <- matrix(data = NA, nrow = 153, ncol = 9)
colnames(res) <- colnames(s)
for(i in seq_along(p[,1])){
  res[i,] <- t(s[i,] - IB %*% A %*% p[i,])
}



# SIGMA0
SIGMA0 <- IB%*%Psi%*%t(IB)+Th

# Get SIGMA ####
SIGMA <- kronecker(RHO, SIGMA0)

plotMat <- function(matrix=matrix){
  plotMat <- t(apply(matrix,2,rev)) # rotate 90 degrees clockwise
  library(lattice)
  print(levelplot(plotMat, ylab="row", xlab="column"))
}
plotMat(RHO)
plotMat(SIGMA0)
plotMat(SIGMA[1:120,1:120])

# Optimization ####
# function to replace free parameters
new.par.list <- list(A = par.list$beta[1:9,10:18], 
                     B = par.list$beta[1:9,1:9],
                     PSI = par.list$psi[1:9,1:9],
                     alpha = 52,
                     a = 53)

MLIST <- list(A = A, B = B, PSI = PSI, TH = TH, alpha = alpha, a = a)

x2MLIST <- function(x, MLIST) {
  A.x <- x[as.vector(new.par.list$A)[as.vector(new.par.list$A)!=0]]
  B.x <- x[as.vector(new.par.list$B)[as.vector(new.par.list$B)!=0]]
  PSI.x <- x[as.vector(new.par.list$PSI)[as.vector(new.par.list$PSI)!=0]]
  alpha <- x[52]
  a <- x[53]
  
  MLIST$A[which(as.vector(new.par.list$A)!=0)]          <- A.x
  MLIST$B[which(as.vector(new.par.list$B)!=0)]          <- B.x
  MLIST$PSI[which(as.vector(new.par.list$PSI)!=0)]      <- PSI.x
  MLIST$alpha[which(as.vector(new.par.list$alpha)!=0)]  <- alpha.x
  MLIST$a[which(as.vector(new.par.list$a)!=0)]          <- a.x
  MLIST
}
# function
# get residuals
get.res <- function (MLIST = NULL, y = y, x = x){
  A <- MLIST$A
  B <- MLIST$B
  x <- p
  y <- s
  res <- matrix(data = NA, nrow = 153, ncol = 9)
  for(i in seq_along(p[,1])){
    res[i,] <- t(s[i,] - IB %*% A %*% p[i,])
  }
  colnames(res) <- colnames(y)
  res
}
head(get.res(MLIST = MLIST, y = s, x = p)) 
# 
get.IB.inv <- function (MLIST = NULL) {
  BETA <- MLIST$B
  nr <- nrow(MLIST$PSI)
  tmp <- -BETA
  tmp[lav_matrix_diag_idx(nr)] <- 1
  IB.inv <- solve(tmp)
}


computeSigmaHat.LISREL <- function (MLIST = NULL) 
{
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  PSI <- MLIST$psi
  THETA <- MLIST$theta
  BETA <- MLIST$beta
  library(lavaan)
  if (is.null(BETA)) {
    LAMBDA..IB.inv <- LAMBDA
  }
  else {
    lavaan:::.internal_get_IB.inv(MLIST = MLIST)
    IB.inv <- get.IB.inv(MLIST = MLIST)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }
  VYx <- tcrossprod(LAMBDA..IB.inv %*% PSI, LAMBDA..IB.inv) + # faster than t(x) %*% y
    THETA
  if (delta && !is.null(MLIST$delta)) {
    DELTA <- diag(MLIST$delta[, 1L], nrow = nvar, ncol = nvar)
    VYx <- DELTA %*% VYx %*% DELTA
  }
  VYx
}

# objective function 'ML'

# Sigma <- matrix(data = c(0.59156973, 0.46014155, 0.65501854, 0.05000000, 0.05664468, 0.04620913, 0.05000000, 0.06125368, 0.04272028, 0.46014155, 0.35791257, 0.50949402, 0.03889157, 0.04406002, 0.03594292, 0.03889157, 0.04764504, 0.03322918, 0.65501854, 0.50949402, 0.72527254, 0.05536275, 0.06272011, 0.05116529, 0.05536275, 0.06782344, 0.04730225, 0.05000000, 0.03889157, 0.05536275, 1.19017633, 0.57889891, 0.47224938, 0.05000000, 0.06125368, 0.04272028, 0.05664468, 0.04406002, 0.06272011, 0.57889891, 1.34672282, 0.53500831, 0.05664468, 0.06939390, 0.04839754, 0.04620913, 0.03594292, 0.05116529, 0.47224938, 0.53500831, 1.07387710, 0.04620913, 0.05660959, 0.03948134, 0.05000000, 0.03889157, 0.05536275, 0.05000000, 0.05664468, 0.04620913, 1.18283419, 0.62172721, 0.43361254, 0.06125368, 0.04764504, 0.06782344, 0.06125368, 0.06939390, 0.05660959, 0.62172721, 1.59155446, 0.53120726, 0.04272028, 0.03322918, 0.04730225, 0.04272028, 0.04839754, 0.03948134, 0.43361254, 0.53120726, 0.96866021),
#                 nrow = 9, ncol = 9)
objective_ML <- function(x, MLIST) {
  MLIST <- x2MLIST(x = x, MLIST = MLIST)
  # compute Sigma.hat
  Sigma <- computeSigmaHat.LISREL(MLIST = MLIST)
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
# get starting values
start.x <- parTable(fit)$start[ parTable(fit)$free > 0 ]

# estimate parameters ML
out.ML  <- nlminb(start = start.x, objective = objective_ML, 
                  MLIST = MLIST, control = list(eval.max = 5000, trace = 1))

x2MLIST(out.ML$par, MLIST)


