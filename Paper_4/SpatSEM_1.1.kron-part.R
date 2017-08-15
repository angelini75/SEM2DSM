rm(list=ls()[])
library(lavaan)

# load lavaan model (from RS server)
setwd("~/big/SEM2DSM1/Paper_4/data")

# load lavaan model (from RS desktop)
setwd("~/Documents/SEM2DSM1/Paper_4/data")

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
#SIGMA.all <- kronecker(SIGMA0, RHO)  # create covariance matrix of z.all
#plotMat(SIGMA.all[1000:1700,1000:1700])
# SIGMA.all2 <- kronM(RHO = RHO,RHO.I = RHO,SIGMA0 = SIGMA0, sp = 1:9, cov = 10:18)
# plotMat(SIGMA.all1[1000:1500,1000:1500]-SIGMA.all2[1000:1500,1000:1500])


# Choleschy decomposition
# L.all = chol(SIGMA.all)
# logdetSIGMA.all = 2*sum(log(diag(L.all)))
# SIGMA.all.inv <- chol2inv(L.all)

p=18
N=153



system.time(c(
  SIGMA0.inv <- chol2inv(chol(SIGMA0)),
   RHO.inv <- chol2inv(chol(RHO)),
   SIGMA.all.inv <- kronecker(SIGMA0.inv, RHO.inv),
   dL.S.R <- append((diag(chol(SIGMA0)))^153, (diag(chol(RHO)))^18),
   logdetSIGMA.all = 2*sum(log(dL.S.R))
))
system.time(c(
  SIGMA.all <- kronecker(SIGMA0, RHO),
  L.all = chol(SIGMA.all),
  logdetSIGMA.all = 2*sum(log(diag(L.all))),
  SIGMA.all.inv <- chol2inv(L.all)
))

# to define the objective function:
objective_ML <- function(x, MLIST = MLIST) {
  MLIST <- x2MLIST(x = x, MLIST = MLIST)
  # compute Sigma.hat
  SIGMA0 <- computeSigmaHat.LISREL(MLIST = MLIST)
  RHO <- get.RHO(MLIST,h)
  if (all(eigen(SIGMA0)$values >0) & (all(eigen(RHO)$values >0))) {
    # L.SIGMA0 <- chol(SIGMA0)
    # L.RHO <- chol(RHO)
    SIGMA0.inv <- chol2inv(chol(SIGMA0))
    RHO.inv <- chol2inv(chol(RHO))
    SIGMA.all.inv <- kronecker(SIGMA0.inv, RHO.inv)
    dL.S.R <- append((diag(chol(SIGMA0)))^153, (diag(chol(RHO)))^18)
    logdetSIGMA.all = 2*sum(log(dL.S.R))
    # SIGMA.all <- kronecker(SIGMA0, RHO)
    # L.all = chol(SIGMA.all)
    # logdetSIGMA.all = 2*sum(log(diag(L.all)))
    # SIGMA.all.inv <- chol2inv(L.all)
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

start.x <- c(lav.est, 0.4, 0.4)
sp.ou <- nlminb(start = start.x, objective = objective_ML, 
                MLIST = MLIST, control = list(iter.max = 500, trace = 1))

#round((start.x - sp.ou$par),4)
MLIST.out <- x2MLIST(sp.ou$par, MLIST)

library(numDeriv)
jacobian <- jacobian(objective_ML, x=sp.ou$par, MLIST=MLIST.out)

# sp.ou2 <- nlminb(start = start.x, objective = objective_ML, 
#                 MLIST = MLIST, control = list(iter.max = 500, trace = 1, 
#                                               rel.tol = 1e-14, x.tol = 1e-12))
rm(list=ls()[])
load("~/Documents/SEM2DSM1/Paper_4/data/9August_SpatSEM_1.1.RData")
################################################################################
# prediction ####
################################################################################
library(gstat)
library(sp)

# function to obtain residuals
# get.res <- function (m = NULL){
#   A <- m$beta[1:9,10:18]
#   B <- m$beta[1:9,1:9]
#   I <- diag(nrow = 9, ncol = 9)
#   IB.inv <- solve(I - B)
#   sp <- z[,1:9]
#   p <- z[,10:18]
#   res <- matrix(data = NA, nrow = 153, ncol = 9)
#   for(i in seq_along(p[,1])){
#     res[i,] <- t(sp[i,] - (IB.inv %*% A %*% p[i,]))
#   }
#   colnames(res) <- colnames(sp)
#   res
# }
## compute the residuals from multivariate linear model
# r <- get.res(m = MLIST.out)
# res <- ks
# res@data[,2:10] <- r
# res@data <- res@data[,2:10]
#var.res <- apply(X = r,FUN =  var, MARGIN = 2)

ks <- as.data.frame(ks)
ks$X <- (ks$X * ST$std.dev[20]) + ST$mean[20]
ks$Y2 <- (ks$Y2 * ST$std.dev[20]) + ST$mean[21]
coordinates(ks) <- ~X+Y2


# spatial model
cv <- gstat(id = "CEC.A", formula = CEC.A ~ 1, data = ks, nmax = 10)
cv <- gstat(cv, id = "CEC.B", formula = CEC.B ~ 1, data = ks, nmax = 10)
cv <- gstat(cv, id = "CEC.C", formula = CEC.C ~ 1, data = ks, nmax = 10)
cv <- gstat(cv, id = "OC.A", formula = OC.A ~ 1, data = ks, nmax = 10)
cv <- gstat(cv, id = "OC.B", formula = OC.B ~ 1, data = ks, nmax = 10)
cv <- gstat(cv, id = "OC.C", formula = OC.C ~ 1, data = ks, nmax = 10)
cv <- gstat(cv, id = "Clay.A", formula = Clay.A ~ 1, data = ks, nmax = 10)
cv <- gstat(cv, id = "Clay.B", formula = Clay.B ~ 1, data = ks, nmax = 10)
cv <- gstat(cv, id = "Clay.C", formula = Clay.C ~ 1, data = ks, nmax = 10)
cv <- gstat(cv, # To fill in the object
            model = vgm(nugget = 99 ,
                        psill= 100,
                        range= MLIST.out$a * ST$std.dev[20], # estimated range times 
                        model = "Exp"), 
            fill.all = T)
# replace nugget and psill # [[[[[[[[[[[[to be continued...]]]]]]]]]]]]]]]]]]]]]
sigma0 <- computeSigmaHat.LISREL(MLIST = MLIST.out)[1:9,1:9]
psi <- MLIST.out$psi[1:9,1:9] # system variance-covariance matrix
#variance <- sigma0-psi
psill <- sigma0[lower.tri(sigma0,T)] * (1-MLIST.out$alpha) - psi[lower.tri(psi,T)]# psill for those free variables 
nugget <- sigma0[lower.tri(sigma0,T)] * MLIST.out$alpha # c0 = psill/alpha

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
png(filename = "~/Dropbox/PhD_Marcos/Paper 4/TeX/Figures/Semivar.png",
width = 3000, height = 3000, res =  250)
plot(cv.var, model=cv.fit, 
     main="",
     xlab = "Distance (m)", 
     ylab = "Semivariance",
     ylim = c(0,1.5),
     scales=list(x = list(alternating = 1), y = list(alternating = 1)),
     par.settings=list(grid.pars=list(fontfamily="serif")))
dev.off()

NAD83.KS.N <- CRS("+init=epsg:2796")
# Assign projection
proj4string(res) <- NAD83.KS.N
mapview(res)

# map of predictions
# libraries
library(raster)
library(maptools)
library(sp)
#library(sf)
library(rgdal)

#define crs
wgs84 <- CRS("+init=epsg:4326")
NAD83.KS.N <- CRS("+init=epsg:2796")
modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# points based on MODIS LST grid 
r <- raster("/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling/LST.mean.tif")
plot(r)
r <- rasterToPoints(r,spatial = TRUE)
proj4string(r) <- modis

# study area WGS84
area <- readOGR("/mnt/L0135974_DATA/UserData/BaseARG/study area/USA/Platte_area_extended.shp")
proj4string(area) <- wgs84
area <- spTransform(area, modis)
area <- as(area,"SpatialPolygons")

# remove points outside the study area
r <- r[!is.na(over(r,area)),]
r <- spTransform(r, NAD83.KS.N)
spplot(r)

# wd
setwd("/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling")
# sdat files (dem covariates) 
files <- list.files(pattern=".sdat$")
header <- gsub(".sdat", "", files)
header <- c("dem", "twi","vdchn") 

# tif files (modis)
files_m <- list.files(pattern=".tif$")
# set names of covariates
header_m <- c("evim", "evisd", "lstm", "lstsd", "ndwi.a", "ndwi.b")

# extract values from files (.sdat)
stack <- list()
for(i in seq_along(files)) {
  r@data[,i] <- NULL
  stack[[i]] <- readGDAL(files[i])
  proj4string(stack[[i]]) <- NAD83.KS.N
  r@data[,i] <- over(r, stack[[i]])[,1]
  stack <- list()
  names(r@data)[length(r@data)] <- header[i]
}  

## extract values from modis files 
stack <- list()
# reproject endo to modis projection
r <- spTransform(r, modis)
for(i in seq_along(files_m)) {
  r@data[,length(r@data)+1] <- NULL
  stack[[i]] <- readGDAL(files_m[i])
  proj4string(stack[[i]]) <- modis # change projection
  r@data[,length(r@data)+1] <- over(r, stack[[i]])[,1]
  stack <- list()
  names(r@data)[length(r@data)] <- header_m[i]
}  
r <- spTransform(r, NAD83.KS.N)
r.df <- as.data.frame(r)
r.df <- r.df[complete.cases(r.df),] 

# transform
r.df$twi <- log10(r.df$twi)
r.df$vdchn <- log10(r.df$vdchn+10)
r.df$ndwi.a <- (r.df$ndwi.a+10)^.3

# standardise #
std <- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] - st[i,2]) / st[i,3]
  }
  y
}
STt <- read.csv("~/Documents/SEM2DSM1/Paper_4/data/STt.ks.csv")
r.st <- std(x = r.df,st = STt[11:21,])
name(pred.st)
