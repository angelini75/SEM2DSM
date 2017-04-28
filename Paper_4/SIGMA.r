library(lattice)

# 9 soil properties
m <- 9
# 7 covariates
q <- 7
# 153 samples
n=153 # n=153 in reality

# Get alpha ####
# alpha = C/(C0 + C) = 0.8/(0.2 + 0.8) 
alpha <- 0.8 
# Get range ####
a <- 1e+05

# Get distances (h) ####
# Load locations of profiles
setwd("~/Documents/SEM2DSM1/Paper_4/data") # change folder
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
# The starting values of Sigma zero (SIGMA0) 
# Matrix of Beta coefficients
#(B <- lavTech(my.fit.lv.ML, "start")$beta[1:9,1:9])
B <- matrix(0, nrow = 9, ncol = 9)
# Identity matrix
I <- diag(nrow = 9, ncol = 9)
# Matrix of Gamma coefficients
#A <- lavTech(my.fit.lv.ML, "start")$beta[1:9,10:18]
A <- matrix(0, nrow = 9, ncol = 9) # 9 covariates x 9 soil properties
# Matrix of Psi coefficients (model error variance-covariance)
#(Psi <- inspect(my.fit.lv.ML, "start")$psi[1:9,1:9])
Psi <- diag(0.05, nrow = 9,ncol = 9)
# Matrix of measurement error (Epsylon)
#Th <- lavTech(my.fit.lv.ML, "start")$theta[1:9,1:9]
Th <- diag(c(0.05,0.05,0.05,0.1,0.1,0.1,0.05,0.05,0.05), 9,9)
IB <- solve(I - B)

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
