library(lavaan)
library(lattice)
library(sp)
library(gstat)
library(pheatmap)

rm(list=ls()[ls()!="my.fit.lv.ML"])
# utils


#plot a heatmap using the calculated correlation matrix
name <- function(x) { as.data.frame(names(x))}
plotMat <- function(matrix=matrix){
  library(lattice)
  a <- levelplot(matrix, ylab="row", xlab="column")#, 
  #              at = seq(from=-1 , to=1, by=0.1)))
  b <- update(a, ylim = rev(extendrange(seq_along(matrix[,1]), f=0.01)))
  b
}
################################################################################
################ HOW LAVAAN ESTIMATES THE COEFICIENTS? #########################
################################################################################
# load lavaan model 
setwd("~/Documents/SEM2DSM1/Paper_4/data")
load("lavaan.model.RData")

# Load data with locations in the same order that is defined in lavaan
names <- rownames(inspect(my.fit.lv.ML, "est")$lambda)[1:18]
ks <- read.csv("ks.csv")[,-1] # standardized data
#ks <- ks[c(-151,-152),]
s <- as.matrix(ks[,names]) 

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

# lavMLIST "start" have the starting values of free parameters
lavMLIST <- lavTech(my.fit.lv.ML, "est")
lavMLIST$a <- a
lavMLIST$alpha <- alpha
# get.res <- function (m = lavMLIST){
#   A <- m$beta[1:9,10:18]
#   B <- m$beta[1:9,1:9]
#   I <- diag(nrow = 9, ncol = 9)
#   IB.inv <- solve(I - B)
#   p <- s[,10:18]
#   v <- s[,1:9]
#   res <- matrix(data = NA, nrow = 153, ncol = 9)
#   for(i in seq_along(p[,1])){
#     res[i,] <- t(v[i,] - (IB.inv %*% A %*% p[i,]))
#   }
#   colnames(res) <- colnames(s)[1:9]
#   res
# }
# res <- get.res(lavMLIST)
# res <- cbind(res, s[,c(12,17)])
# res <- as.data.frame(res)
# coordinates(res) <- ~X+Y

cv <- gstat(id = "CEC.A", formula = eval(parse(text = paste0(names[1]," ~ 1"))), 
            data = ks, nmax = 10)
for( i in 2:18){
  cv <- gstat(cv, id = names[i], formula = eval(parse(text = paste0(names[i]," ~ 1"))),
              data = ks, nmax = 10)
}
cv <- gstat(cv, 
            model = vgm(nugget = 0.30,
                        psill= 1.4,
                        range=0.8,
                        model = "Exp"), 
            fill.all = T)
cv
cut <- 1
wid <- 0.05
cv.var<- variogram(object = cv, cutoff = cut, width = wid, covariogram = TRUE) 


#cv.fit<-fit.lmc(v = cv.var,g =  cv, fit.lmc = F, fit.ranges = F) 

# png(filename = "~/Dropbox/PhD Marcos/Paper 4/Figures/Fig2.png", 
# width = 3000, height = 3000, res =  250)
plot(cv.var, #model=cv.fit, 
     main="Variograms and cross-variograms of standardized residuals",
     xlab = "Distance / m", 
     ylab = "Semivariance", cex=0.2,
     scales=list(x = list(alternating = 1), y = list(alternating = 1)),
     par.settings=list(grid.pars=list(fontfamily="serif")))


cv.var <- as.data.frame(cv.var)
cv.mat <- matrix(cv.var$gamma, 21)
cv.mat <- rbind(cv.mat[21,],cv.mat[1:20,])
colnames(cv.mat) <- unique(cv.var$id)
cv.mat <- as.data.frame(cv.mat)
cv.mat$dist.l <- seq(from = 0, to = (cut), by = wid)
cv.mat$dist.u <- seq(from = wid, to = cut+0.05, by = wid)
names(cv.mat)[which(names(cv.mat) %in% names)] <- 
  paste(names,names, sep = ".")
#cv.mat <- cv.mat[c(1,46,47)]
cv.mat$dist.u[which(cv.mat$dist.u == max(cv.mat$dist.u))] <- 5

sample <- matrix(1, 18,18)
S <- kronecker(h,sample)
plotMat(S[1:360,1:360])
ev <- matrix(NA, 153, 153)
for(k in 1:153){ # location 1, 2, 3, ... , k.
  for(l in 1:153){ # location 1, 2, 3, ... , l.
    dist <- h[l,k] # distance between the pair of points k and l.
    if (dist != 0) { # for k != l
      s <- matrix(NA, 18,18) # create a matrix of (q+p) x (q+p)
      # dist have to match a "d" lag distance at matrix cv.mat 
      # (d is the same for the whole s matrix)
      d <- which(dist >= cv.mat$dist.l & dist < cv.mat$dist.u) 
      for(j in seq_along(names)){ # variable j within the 18x18 matrix
        for(i in seq_along(names)){ # variable i within the 18x18 matrix
          if (j > i) { # as it is a symmetric matrix, s[i,j] = s[j,i]
            var.var <- paste(names[i], names[j], sep = ".")
          } else {
            var.var <- paste(names[j], names[i], sep = ".")
          }
          # sample variance for var.var pair and d lag distance 
          gamma <- cv.mat[d,var.var] 
          s[i,j] <- gamma # at location i,j of s matrix
        }
      }
    } else { # if dist = 0, s = SIGMA0
      s <- my.fit.lv.ML@SampleStats@cov[[1]]
    } # s is inserted in lo
    S[((l*18-17):(l*18)),((k*18-17):(k*18))] <- s
  }
}
plotMat(s)
plotMat(S[1:360,1:360])

get.S <- function(h, cv.mat){
  m <- h
  for (i in seq_along(h){
    dist <- which(h[2] > cv.mat$dist.l & h[2] < cv.mat$dist.u)
    m <- cv.mat[dist,1]
  }
  m
}

#






# adapted from https://groups.google.com/forum/#!searchin/lavaan/lavaan$20function$20github%7Csort:relevance/lavaan/0Hitqi3k_do/pYfyoocABwAJ

# This are the estimates of our model
lavTech(my.fit.lv.ML, "est")

# lavMLIST "start" have the starting values of free parameters
lavMLIST <- lavTech(my.fit.lv.ML, "start")
lavMLIST$a <- a
lavMLIST$alpha <- alpha
# This is the sample covariance matrix called S
#S <- my.fit.lv.ML@SampleStats@cov[[1]]
#S <- kronecker(RHO, S)
# inverse and log(det(S))
cS <- chol(S)
S.inv <- chol2inv(cS)
diag.cS <- diag(cS)
S.logdet <- sum(log(diag.cS * diag.cS))

# x2MLIST (renamed as x2lavMLIST) fill in new values (from x) in the 
# model matrices (lavMLIST)
# lavaan:::MLISTX2MLIST

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
get.SIGMA0 <- function (lavMLIST = NULL)#, delta = TRUE) 
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
  # if (delta && !is.null(lavMLIST$delta)) {
  #   DELTA <- diag(lavMLIST$delta[, 1L], nrow = nvar, ncol = nvar)
  #   VYx <- DELTA %*% VYx %*% DELTA
  # }
  VYx
}

get.RHO <- function(lavMLIST = NULL) {
  a <- lavMLIST$a
  n <- nrow(h)
  alpha <- lavMLIST$alpha
  RHO <- matrix(rep(NA,n^2), nrow = n)
  for(i in seq_along(RHO)) {
    RHO[i] <- alpha * exp(-h[i]/a)
  }
  diag(RHO) <- 1
  RHO
}
plotMat(get.RHO(lavMLIST = lavMLIST)[1:27,1:27])

get.SIGMA <- function(SIGMA0 = NULL, RHO = NULL){
  V <- kronecker(SIGMA0, RHO, "*")
  V
}

SIGMA0 <- get.SIGMA0(lavMLIST)
RHO <- get.RHO(lavMLIST)
plotMat(get.SIGMA(SIGMA0, RHO)[1:90,1:90])
# objective function 'ML'
