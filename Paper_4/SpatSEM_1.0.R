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


################################################################################
# prediction ####
################################################################################
library(gstat)
library(sp)

# function to obtain residuals
get.res <- function (m = NULL){
  A <- m$beta[1:9,10:18]
  B <- m$beta[1:9,1:9]
  I <- diag(nrow = 9, ncol = 9)
  IB.inv <- solve(I - B)
  sp <- z[,1:9]
  p <- z[,10:18]
  res <- matrix(data = NA, nrow = 153, ncol = 9)
  for(i in seq_along(p[,1])){
    res[i,] <- t(sp[i,] - (IB.inv %*% A %*% p[i,]))
  }
  colnames(res) <- colnames(sp)
  res
}
## compute the residuals from multivariate linear model
r <- get.res(m = MLIST.out)
res <- ks
res@data[,2:10] <- r
res@data <- res@data[,2:10]
#var.res <- apply(X = r,FUN =  var, MARGIN = 2)

res <- as.data.frame(res)
res$X <- (res$X * ST$std.dev[20]) + ST$mean[20]
res$Y2 <- (res$Y2 * ST$std.dev[20]) + ST$mean[21]
coordinates(res) <- ~X+Y2


# spatial model
cv <- gstat(id = "CEC.A", formula = CEC.A ~ 1, data = res, nmax = 10)
cv <- gstat(cv, id = "CEC.B", formula = CEC.B ~ 1, data = res, nmax = 10)
cv <- gstat(cv, id = "CEC.C", formula = CEC.C ~ 1, data = res, nmax = 10)
cv <- gstat(cv, id = "OC.A", formula = OC.A ~ 1, data = res, nmax = 10)
cv <- gstat(cv, id = "OC.B", formula = OC.B ~ 1, data = res, nmax = 10)
cv <- gstat(cv, id = "OC.C", formula = OC.C ~ 1, data = res, nmax = 10)
cv <- gstat(cv, id = "Clay.A", formula = Clay.A ~ 1, data = res, nmax = 10)
cv <- gstat(cv, id = "Clay.B", formula = Clay.B ~ 1, data = res, nmax = 10)
cv <- gstat(cv, id = "Clay.C", formula = Clay.C ~ 1, data = res, nmax = 10)
cv <- gstat(cv, # To fill in the object
            model = vgm(nugget = 0 ,
                        psill= 0,
                        range= MLIST.out$a * ST$std.dev[20], # estimated range times 
                        model = "Exp"), 
            fill.all = T)
# replace nugget and psill # [[[[[[[[[[[[to be continued...]]]]]]]]]]]]]]]]]]]]]
sigma0 <- computeSigmaHat.LISREL(MLIST = MLIST.out)[1:9,1:9]
psi <- MLIST.out$psi[1:9,1:9] # system variance-covariance matrix
variance <- sigma0-psi
psill <- variance[lower.tri(variance,T)] # psill for those free variables 
nugget <- variance[lower.tri(variance,T)] * MLIST.out$alpha # c0 = psill/alpha

# replace nugget and psill values in cv object
for(i in 1:45){
  cv$model[[i]][1,"psill"] <- nugget[i]
  cv$model[[i]][2,"psill"] <- psill[i]
}

#
cv.var<- variogram(object = cv, cutoff = 150000) 
cv.fit<-fit.lmc(v = cv.var,g =  cv, fit.lmc = F, fit.ranges = F) 

# 
cv.fit$model <- cv$model
# png(filename = "~/Dropbox/PhD Marcos/Paper 4/Figures/Fig2.png", 
# width = 3000, height = 3000, res =  250)
plot(cv.var, model=cv.fit, 
     main="Variograms and cross-variograms of standardized residuals",
     xlab = "Distance / m", 
     ylab = "Semivariance",
     scales=list(x = list(alternating = 1), y = list(alternating = 1)),
     par.settings=list(grid.pars=list(fontfamily="serif")))
# dev.off()
NAD83.KS.N <- CRS("+init=epsg:2796")
# Assign projection
proj4string(res) <- NAD83.KS.N
mapview(res)

# > sigma0
# [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]      [,9]
# [1,] 1.3341429 0.8517236 0.5566962 0.4038747 0.3496910 0.2953360 1.1694629 0.7503792 0.6028779
# [2,] 0.8517236 1.2546326 0.8225599 0.3285202 0.4073502 0.1803358 0.7540358 1.0353120 0.7111746
# [3,] 0.5566962 0.8225599 1.1344237 0.1867062 0.1378134 0.1224321 0.4932839 0.6215340 0.8643812
# [4,] 0.4038747 0.3285202 0.1867062 1.2046502 0.4991251 0.3728067 0.2123006 0.3631923 0.2477451
# [5,] 0.3496910 0.4073502 0.1378134 0.4991251 1.5236514 1.1729970 0.2728593 0.4183592 0.2292178
# [6,] 0.2953360 0.1803358 0.1224321 0.3728067 1.1729970 1.4480942 0.2424055 0.1817066 0.2036348
# [7,] 1.1694629 0.7540358 0.4932839 0.2123006 0.2728593 0.2424055 1.2677797 0.8239197 0.6178310
# [8,] 0.7503792 1.0353120 0.6215340 0.3631923 0.4183592 0.1817066 0.8239197 1.1848127 0.7863506
# [9,] 0.6028779 0.7111746 0.8643812 0.2477451 0.2292178 0.2036348 0.6178310 0.7863506 1.1113939
# 
# > psi
# [,1]      [,2]       [,3]      [,4]     [,5]       [,6]      [,7]        [,8]      [,9]
# [1,]  0.16722111 0.1607945 0.08297455 0.0000000 0.000000  0.0000000 0.0000000 -0.07354172 0.0000000
# [2,]  0.16079454 0.2933525 0.25223058 0.0000000 0.000000  0.0000000 0.0000000  0.00000000 0.0000000
# [3,]  0.08297455 0.2522306 0.39957203 0.0000000 0.000000  0.0000000 0.0000000  0.00000000 0.0000000
# [4,]  0.00000000 0.0000000 0.00000000 0.9748911 0.000000  0.0000000 0.0000000  0.00000000 0.0000000
# [5,]  0.00000000 0.0000000 0.00000000 0.0000000 0.971422  0.0000000 0.0000000  0.00000000 0.0000000
# [6,]  0.00000000 0.0000000 0.00000000 0.0000000 0.000000  0.3875347 0.0000000 -0.18995953 0.0000000
# [7,]  0.00000000 0.0000000 0.00000000 0.0000000 0.000000  0.0000000 0.7714742  0.00000000 0.0000000
# [8,] -0.07354172 0.0000000 0.00000000 0.0000000 0.000000 -0.1899595 0.0000000  0.35693807 0.0000000
# [9,]  0.00000000 0.0000000 0.00000000 0.0000000 0.000000  0.0000000 0.0000000  0.00000000 0.8732465
