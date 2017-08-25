# LOOCV 

rm(list=ls()[])
library(lavaan)
library(doParallel)
library(sp)

# load lavaan model (from RS server)
setwd("~/big/SEM2DSM1/Paper_4/data")

# load lavaan model (from RS desktop)
setwd("~/Documents/SEM2DSM1/Paper_4/data")
#setwd("C:/Users/quics/Marcos/SEM2DSM1/Paper_4/data")
load("SpatSEM_1.2.RData")
# ks <- read.csv("ks.csv")
# ks <- ks[,colnames(s)]
# s <- as.matrix(ks[,-17])
# 
# # lavaan model
# fit <- my.fit.lv.ML
# par.list <- inspect(fit)
# # estimated parameters with lavaan
# MLIST <- lavTech(fit, "est")
# 
# # compute SIGMA0 (18x18)
# SIGMA0 <- computeSigmaHat.LISREL(MLIST)
# #plotMat(SIGMA0)
# 
# # initialise matrix with standardised observations
# # rows are locations, columns variables
# z <- as.matrix(s)
# 
###############################################################
########### FUNCTIONS              ###########################
#############################################################
# standardize variables 
std <- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] - st[i,2]) / st[i,3]
  }
  y
}
# Function to include free parameters into the matrix structure MLIST
x2MLIST <- function(x, MLIST) {
  lambda.x <- x[as.vector(par.list$lambda)[as.vector(par.list$lambda)!=0]]
  theta.x <- x[as.vector(par.list$theta)[as.vector(par.list$theta)!=0]]
  psi.x <- x[as.vector(par.list$psi)[as.vector(par.list$psi)!=0]]
  beta.x <- x[as.vector(par.list$beta)[as.vector(par.list$beta)!=0]]
  alpha.x <- x[free+1]
  a.x <- x[free+2]
  
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

# Objective function:
objective_ML <- function(x, MLIST = MLIST) {
  MLIST <- x2MLIST(x = x, MLIST = MLIST)
  # compute Sigma.hat
  SIGMA0 <- computeSigmaHat.LISREL(MLIST = MLIST)
  RHO <- get.RHO(MLIST,h)
  if (all(eigen(SIGMA0)$values >0) & (all(eigen(RHO)$values >0))) {
    
    SIGMA0.inv <- chol2inv(chol(SIGMA0))
    RHO.inv <- chol2inv(chol(RHO))
    SIGMA.all.inv <- kronecker(SIGMA0.inv, RHO.inv)
    dL.S.R <- append((diag(chol(SIGMA0)))^N, (diag(chol(RHO)))^p)
    logdetSIGMA.all = 2*sum(log(dL.S.R))
    
    objective <- -1 * (-1/2*p*N*log(2*pi) - 1/2*logdetSIGMA.all - 
                         1/2 * crossprod(z.all, SIGMA.all.inv) %*%z.all)
    cat("objective = ", objective, "\n")
    objective
  } else {
    objective <- Inf
    objective
  }
}

# Change coordinates (not for LOOCV)
change.coords <- function(x = samples){
  duplos <- unique(zerodist(x, zero=0.0)[,2])
  coords <- x@coords[duplos,]
  coords.v <- as.vector(as.matrix(coords))
  coords[coords>0] <- coords[coords>0] + 
    rnorm(length(coords[coords>0]), mean = 50, sd = 30)
  coords[coords<0] <- coords[coords<0] +
    rnorm(length(coords[coords<0]), mean = -50, sd = 30)
  x@coords[duplos,] <- coords
  new.coords <- as.data.frame(x@coords)
  x <- as.data.frame(x)
  sp::coordinates(x) <- ~X+Y2
  x
}

# Function to get RHO0 including one prediction location
get.RHO0 <- function(MLIST = NULL, h0 = NULL) {
  a <- MLIST$a
  alpha <- MLIST$alpha
  n <- nrow(h0)
  RHO0 <- matrix(rep(NA,n^2), nrow = nrow(h0), ncol = ncol(h0))
  for(i in seq_along(RHO0)) {
    RHO0[i] <- (1-alpha) * exp(-h0[i]/a)
  }
  #diag(RHO) <- 1
  RHO0
}

## function to get prediction at new location from trend model
# it includes parallel processing
get.pred <- function (MLIST = NULL, covar = covar.st, xy){
  var.names <- c("CEC.Ar","CEC.Br","CEC.Cr","OC.Ar","OC.Br","OC.Cr",
                 "Clay.Ar","Clay.Br","Clay.Cr","dem","vdchn","X",
                 "lstm","evisd","ndwi.b","twi","ndwi.a")
  m <- MLIST
  A <- m$beta[1:9,10:p]
  B <- m$beta[1:9,1:9]
  I <- diag(nrow = 9, ncol = 9)
  IB.inv <- solve(I - B)
  k <- covar[,var.names[10:p]]
  p <- as.vector(as.matrix(k[i,]))
  pred <- t(IB.inv %*% A %*% p) 
  colnames(pred) <- var.names[1:9]
  cbind(pred,xy)
}
# function to get residuals
get.res <- function (m = NULL, z = cal){
  A <- m$beta[1:9,10:p]
  B <- m$beta[1:9,1:9]
  I <- diag(nrow = 9, ncol = 9)
  IB.inv <- solve(I - B)
  sp <- z[,1:9]
  p <- z[,10:p]
  res <- matrix(data = NA, nrow = N, ncol = 9)
  for(i in seq_along(p[,1])){
    res[i,] <- t(sp[i,] - (IB.inv %*% A %*% p[i,]))
  }
  colnames(res) <- colnames(sp)
  res
}

#####################################################################
################ Parameters for LOOCV ##############################
###################################################################

ks <- read.csv("ks.csv")[,c(colnames(s),"Y")] # standardized data
ST <- read.csv("STt.ks-0.3.csv")
rownames(ST) <- ST$X

N = nrow(ks)-1
p = ncol(s)
MLIST <- MLIST
MLIST$alpha <- 0.4
MLIST$a <- 0.4
free <-  length(lavaan:::lav_model_get_parameters(lavmodel = fit@Model))
lav.est <- parTable(fit)$est[parTable(fit)$free > 0]
start.x <- c(lav.est, 0.4, 0.4) #lavaan parameters + alpha + a
predicted <- ks[,1:9]
predicted[,] <- NA
k.residual <- predicted
var.names <- c("CEC.Ar","CEC.Br","CEC.Cr","OC.Ar","OC.Br","OC.Cr",
               "Clay.Ar","Clay.Br","Clay.Cr","dem","vdchn","X",
               "lstm","evisd","ndwi.b","twi","ndwi.a")

#################################################################
############### Foreach loop ###################################
###############################################################
doParallel::registerDoParallel(cores = 12)
system.time({
foreach(i = icount(nrow(12)), .combine = cbind,
        .packages = c("sp")) %dopar% {
  # calibration dataset
  cal <- ks[-i,]
  rownames(cal) <- 1:N
  # get h
  xy <- cal[, c("X","Y")]
  xy[,"Y"] <- xy[,"Y"] * ST$std.dev[21] / ST$std.dev[20]
  coordinates(xy) <- ~X+Y
  h <- sp::spDists(xy)
  RHO <- get.RHO(MLIST = MLIST, h = h)
  cal <- as.matrix(cal[,colnames(s)])
  z.all <- as.vector(cal)
  out <- nlminb(start = start.x, objective = objective_ML, 
                 MLIST = MLIST, control = list(iter.max = 2))
  MLIST.obs <- x2MLIST(out$par, MLIST)
  RHO <- get.RHO(MLIST = MLIST.obs, h = h)
  #
  covar.st <- ks[i,c("dem","vdchn","X", "lstm","evisd","ndwi.b","twi","ndwi.a")]
  ll <- ks[i,c("X","Y")]
  ll[,"Y"] <- ll[,"Y"] * ST$std.dev[21] / ST$std.dev[20]
  SIGMA.yy <- computeSigmaHat.LISREL(MLIST = MLIST.obs)[1:9,1:9] # qxq
  SIGMA.xx <- kronecker(SIGMA.yy, RHO)   # qN x qN
  pred.lm <- get.pred(MLIST = MLIST.obs, covar = covar.st, xy = ll)  
  res <- get.res(m = fit@Model@GLIST)
  y.all <- as.vector(res) # vector of residuals
  xy <- as.data.frame(xy)
  xy.ll <- rbind(xy,ll)
  coordinates(xy.ll) <- ~X+Y
  h.all <- sp::spDists(xy.ll)
  h0 <- matrix(h.all[1:N,N+1], ncol = 1, nrow = N)
  RHO0 <- get.RHO0(MLIST.obs, h0 = h0) # note that h0 is Nx1 
  SIGMA.xy <- kronecker(SIGMA.yy, RHO0) # qN x q
  k.res <- t(crossprod(SIGMA.xy, chol2inv(chol(SIGMA.xx))) %*% y.all)
  colnames(k.res) <- colnames(res)
  # total prediction
  predicted[i,] <- pred.lm + k.res
  k.residual[i,] <- k.res
  # Computing prediction variance
  PSI <- MLIST.obs$psi[1:9,1:9]
  B <- MLIST.obs$beta[1:9,1:9]
  I <- diag(nrow = 9, ncol = 9)
  IB.inv <- solve(I - B)
  theta <- MLIST.obs$theta[1:9,1:9]
  var.SIGMA.yy <- PSI + IB.inv %*% tcrossprod(theta, IB.inv)
  var.SIGMA.xx <- kronecker(var.SIGMA.yy, RHO)
  var.SIGMA.xy <- kronecker(var.SIGMA.yy, RHO0)
  var.zeta <- var.SIGMA.yy - 
                     crossprod(var.SIGMA.xy, 
                               chol2inv(chol(var.SIGMA.xx))) %*% var.SIGMA.xy
  variance[i,] <- diag(var.zeta)
  as.vector(var.zeta)
}
})

(rowMeans(var.zeta)^0.5) * STt$std.dev[c(5:7,8:10,2:4)]
doParallel::stopImplicitCluster()

pred <- pred.lm + k.res
sds <- matrix(rep(STt$std.dev[c(5:7,8:10,2:4)], each=1), ncol=9)
means <- matrix(rep(STt$mean[c(5:7,8:10,2:4)], each=1), ncol=9)
summary((pred[rownames(loc),1:9] + k.res) *  sds + means)
