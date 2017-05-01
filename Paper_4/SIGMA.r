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
get.RHO <- function(range = a, 
                    samples = n, 
                    dist = h, 
                    ratio = alpha ) {
  a <- range
  n <- samples
  h <- dist
  alpha <- ratio
  RHO <- matrix(rep(NA,n^2), nrow = n)
  for(i in seq_along(RHO)) {
    RHO[i] <- alpha * exp(-h[i]/a)
  }
  diag(RHO) <- 1
  RHO
}
RHO <- get.RHO()

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

MLIST <- list(A = A, B = B, PSI = PSI, TH = TH, alpha = alpha, a = a, h = h)

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
# functions ####
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

# IB inv 
get.IB.inv <- function (MLIST = NULL) {
  BETA <- MLIST$B
  nr <- nrow(MLIST$PSI)
  tmp <- -BETA
  tmp[lav_matrix_diag_idx(nr)] <- 1
  IB.inv <- solve(tmp)
}

# 
get.IB.PSI.IBinv <- function (MLIST = NULL) 
{
  PSI <- MLIST$PSI
  THETA <- MLIST$TH
  #BETA <- MLIST$B
  IB.inv <- get.IB.inv(MLIST = MLIST)
  IB.PSI.IBinv <- tcrossprod(IB.inv %*% PSI, IB.inv) + # faster than t(x) %*% y
    THETA
  IB.PSI.IBinv
}

# get RHO
get.RHO <- function(MLIST = NULL, 
                    samples = n) {
  a <- MLIST$a
  n <- nrow(MLIST$h)
  h <- MLIST$h
  alpha <- MLIST$alpha
  RHO <- matrix(rep(NA,n^2), nrow = n)
  for(i in seq_along(RHO)) {
    RHO[i] <- alpha * exp(-h[i]/a)
  }
  diag(RHO) <- 1
  RHO
}
get.RHO(MLIST = MLIST)

# objective function 'ML'

objective_ML <- function(x, MLIST) {
  MLIST <- x2MLIST(x = x, MLIST = MLIST)
  # compute Sigma.hat
  SIGMA0 <- get.IB.PSI.IBinv(MLIST = MLIST)
  RHO <- get.RHO(MLIST = MLIST)
  SIGMA <- kronecker(RHO, SIGMA0)
  if (all(eigen(SIGMA)$values >0)) {
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


