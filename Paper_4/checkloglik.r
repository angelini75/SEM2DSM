rm(list=ls()[])
library(lavaan)

# load lavaan model 
 setwd("~/big/SEM2DSM1/Paper_4/data")
setwd("~/Documents/SEM2DSM1/Paper_4/data")
load("env.for.gerard.RData")
load("lavaan.model.RData")

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
# get the correct h
library(sp)
s <- as.data.frame(s)
coordinates(s) <- ~X+Y2
h <- spDists(s)
s <- as.data.frame(s[,1:18])
s <- as.matrix(s)
#
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
    RHO[i] <- (1-alpha) * exp(-h[i]/a)
  }
  diag(RHO) <- 1
  RHO
}




# KronM function: adaptation of Kronecker function 
#source("~/big/SEM2DSM1/Paper_4/kronM.R")



# next our approach with alpha != 0
MLIST$alpha <- alpha
MLIST$a <- a

z.all <- as.vector(z)  # compile to one big vector
RHO <- get.RHO(MLIST,h)
#SIGMA.all <- kronecker(RHO, SIGMA0)  # create covariance matrix of z.all
SIGMA.all <- kronM(RHO = RHO,RHO.I = diag(153),SIGMA0 = SIGMA0, sp = 1:9)
plotMat(SIGMA.all[1000:1700,1000:1700])

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
  if (all(eigen(SIGMA0)$values >0) & (all(eigen(RHO)$values >0))) {
  SIGMA.all <- kronM(RHO = RHO,RHO.I = diag(153),SIGMA0 = SIGMA0, sp = 1:9, cov = 10:18)
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

# ranges <- c(0.01, 0.02,
#             0.03, 0.04,
#             0.05, 0.06, 
#             0.08, 0.1,
#             0.15, 0.20,
#             0.27, 0.35,
#             0.47, 0.52)
likhood <- numeric()
output <- list()
output[1:14] <- NA
#for(i in seq_along(ranges)){
#  MLIST$a <- ranges[i]
  start.x <- c(lav.est, alpha)
  # optimizer of objective funtion 
  l  <- try(nlminb(start = start.x, objective = objective_ML, 
                    MLIST = MLIST, control = list(iter.max = 200, trace = 1)))
  if ('try-error' %in% class(l)) (output[[i]]=NA & next)
  else output[[i]] <- l
  likhood <- append(likhood, l$objective)
}

start.x <- c(lav.est, 0.7, 0.5)
sp.ou <- nlminb(start = start.x, objective = objective_ML, 
       MLIST = MLIST, control = list(iter.max = 150, trace = 1))

round((start.x - sp.out$par),4)
MLIST.sp <- x2MLIST(sp.out$par, MLIST)

################################################################################
# prediction ####
################################################################################
library(gstat)
# load estimates
load("checkloglik.RData")
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
## compute the residuals from multivariate linear model
r <- get.res(m = MLIST.sp)
res <- ks
res@data[,2:10] <- r
res@data <- res@data[,2:10]
#var.res <- apply(X = r,FUN =  var, MARGIN = 2)

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
                        range= MLIST.sp$a, # estimated range
                        model = "Exp"), 
            fill.all = T)
# replace nugget and psill
psi <- MLIST.sp$psi[1:9,1:9] # system variance-covariance matrix
psill <- psi[lower.tri(psi,T)] # psill for those free variables 
nugget <- psi[lower.tri(psi,T)]/2.077891 # c0 = psill/2.077891

# replace nugget and psill values in cv object
for(i in 1:45){
  cv$model[[i]][1,"psill"] <- nugget[i]
  cv$model[[i]][2,"psill"] <- psill[i]
}

#
cv.var<- variogram(object = cv, cutoff = 0.7) 
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

# location standarized points
xy <- read.csv( "xy.csv")[,-1]
for(i in names(xy)){
  xy[,i] <- (xy[,i] - ST$mean[ST$X == i])/ST$std.dev[ST$X == i]
}
xy$Y2 <- xy$Y * ST$std.dev[ST$X=="Y"]/ST$std.dev[ST$X=="X"]
coordinates(xy) <- ~X+Y2

# prediction function
get.pred <- function (m = NULL, data = s, newdata){
  A <- m$beta[1:9,10:18]
  B <- m$beta[1:9,1:9]
  I <- diag(nrow = 9, ncol = 9)
  IB.inv <- solve(I - B)
  #sp <- s[,1:9]
  p <- newdata[,colnames(s)[10:18]]
  pred <- matrix(data = NA, nrow = nrow(newdata), ncol = 9)
  for(i in seq_along(p[,1])){
    return(i)
    pred[i,] <- t(IB.inv %*% A %*% p[i,])
  }
  colnames(pred) <- colnames(s)[1:9]
  pred
}

# predictions from multivariate linear model
pred.lm <- get.pred(m= MLIST.sp, data = s, newdata = as(xy, "data.frame"))

## krige
data <- ks[,-1]
data@data[1:9] <- res@data
gstat:::predict.gstat(formula = Clay.C ~ 1, data = data, newdata = xy, cv)

#
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y
m <- vgm(.59, "Sph", 874, .04)
# ordinary kriging:
x <- krige(log(zinc)~1, meuse, meuse.grid, model = m)
spplot(x["var1.pred"], main = "ordinary kriging predictions")
spplot(x["var1.var"],  main = "ordinary kriging variance")
# simple kriging:
x <- krige(log(zinc)~1, meuse, meuse.grid, model = m, beta = 5.9)
# residual variogram:
m <- vgm(.4, "Sph", 954, .06)
# universal block kriging:
x <- krige(log(zinc)~x+y, meuse, meuse.grid, model = m, block = c(40,40))
spplot(x["var1.pred"], main = "universal kriging predictions")

v = function(x, y) { exp(-spDists(coordinates(x),coordinates(y))) }


# krige0, using user-defined covariance function and multiple responses in y:
# exponential variogram with range 500, defined as covariance function:
v = function(x, y = x) { exp(-spDists(coordinates(x),coordinates(y))/500) }
# krige two variables in a single pass (using 1 covariance model):
y = cbind(meuse$zinc,meuse$copper,meuse$lead,meuse$cadmium)
x <- krige0(formula = zinc~1,data =  meuse,newdata =  meuse.grid, v, y = y)
meuse.grid$zinc = x[,1]
spplot(meuse.grid["zinc"], main = "zinc")
meuse.grid$copper = x[,2]
spplot(meuse.grid["copper"], main = "copper")





























# accuracy
MLIST.sp <- x2MLIST(sp.out$par, MLIST)

get.pred.sp <- function (m = NULL){
  A <- m$beta[1:9,10:18]
  B <- m$beta[1:9,1:9]
  RHO <- get.RHO(MLIST = m, h)
  I <- diag(nrow = 9, ncol = 9)
  IB.inv <- solve(I - MLIST.sp$beta[1:9,1:9])
  sp <- s[,1:9]
  p <- s[,10:18]
  pred <- matrix(data = NA, nrow = 153, ncol = 9)
  for(i in seq_along(p[,1])){
    pred[i,] <- IB.inv %*% A %*% p[i,]
  }
  colnames(res) <- colnames(sp)
  hist((t(RHO) %*% (sp - pred))[,4])
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
