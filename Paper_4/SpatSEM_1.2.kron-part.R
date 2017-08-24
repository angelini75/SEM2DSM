rm(list=ls()[])
library(lavaan)
library(doParallel)
library(sp)

# load lavaan model (from RS server)
setwd("~/big/SEM2DSM1/Paper_4/data")

# load lavaan model (from RS desktop)
setwd("~/Documents/SEM2DSM1/Paper_4/data")
setwd("C:/Users/quics/Marcos/SEM2DSM1/Paper_4/data")
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
# ################################################################################
ks <- read.csv("ks.csv")[,-1] # standardized data
ST <- read.csv("STt.ks-0.3.csv")
ks$Y2 <- ks$Y * ST$std.dev[21] / ST$std.dev[20]

coordinates(ks) <- ~X+Y2
# h = n x n matrix of distances between samples 
h <- sp::spDists(ks)


# Define the parameters alpha and range, and estimate distance matrix:
# Sill to nugget ratio (alpha) ####
# alpha = C/(C0 + C) = 0.8/(0.2 + 0.8) 
alpha <- 0.5
# Range (a) ####
a <- 0.7 #e+05
# number of variables
p=17
# number of samples
N=147
# number of free parameters of lavaan model
free <-  length(lavaan:::lav_model_get_parameters(lavmodel = fit@Model))

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

# KronM function: adaptation of Kronecker function [NOT NECESSARY]
#source("~/big/SEM2DSM1/Paper_4/kronM.R")

# next our approach with alpha != 0
MLIST$alpha <- alpha
MLIST$a <- a

z.all <- as.vector(s)

# to define the objective function:
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
                         1/2*t(z.all)%*%SIGMA.all.inv%*%z.all)
    cat("objective = ", objective, "\n")
    objective
  } else {
    objective <- Inf
    objective
  }
}



# get the 51 starting values to feed x
lav.est <- parTable(fit)$est[parTable(fit)$free > 0]
lavaan:::lav_model_get_parameters(lavmodel = fit@Model)
# try with fixted a and alpha
start.x <- c(lav.est, 0.4, 0.4)

#### BOOTSTRAPPING #####
#registerDoParallel(cores = 46)
trials = 1
#ptime <- system.time({

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

# x <- foreach(icount(trials), .combine=rbind) %dopar% {
#   
# samples <- as.data.ks#[sample(x = 1:nrow(ks), size = nrow(ks), replace = TRUE),]
# rownames(samples) <- 1:147
# samples$Y2 <- samples$Y * ST$std.dev[21] / ST$std.dev[20]
# sp::coordinates(samples) <- ~X+Y2
# samples.ch <- change.coords(x = samples)
# 
# h <- sp::spDists(samples)#.ch)
# #plotMat(h)
# RHO <- get.RHO(MLIST,h)
# #plotMat(RHO)
# z <- as.matrix(as.data.frame(samples)[,colnames(s)])
# z.all <- as.vector(z)
# # out <-  try(nlminb(start = start.x, objective = objective_ML, 
# #                    MLIST = MLIST, control = list(trace = 1)), silent = TRUE)
# #   if(inherits(out, "try-error")) {NA} else {
# #     out$par
#   }
# }
out <-  nlminb(start = start.x, objective = objective_ML, 
               MLIST = MLIST, control = list(trace = 1))

MLIST.out <- x2MLIST(out$par, MLIST)



######### Matrices for prediction 
# for x = zeta_obs (s_i) and y = zeta_obs (s_0)
MLIST.obs <- MLIST.out
# MLIST.obs$lambda <- MLIST.obs$lambda[1:9,1:9]
# MLIST.obs$theta <- MLIST.obs$theta[1:9,1:9]
# MLIST.obs$psi <- MLIST.obs$psi[1:9,1:9]
# MLIST.obs$beta <- MLIST.obs$beta[1:9,1:9]

# load and standardize prediction locations 
std <- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] - st[i,2]) / st[i,3]
  }
  y
}
# standardize
covar <- read.csv("~/big/SEM2DSM1/Paper_4/data/xy.csv")[,-1]
STt <- read.csv("~/big/SEM2DSM1/Paper_4/data/STt.ks-0.3.csv")
covar.st <- std(x = covar[,as.character(STt[c(11:13,15,16,18:21),1])],
                st = STt[c(11:13,15,16,18:21),])
covar.st <- covar.st[complete.cases(covar.st),]

# Function to get RHO including one prediction location
get.RHO0 <- function(MLIST = NULL, h = h) {
  a <- MLIST$a
  n <- nrow(h0)
  alpha <- MLIST$alpha
  RHO <- matrix(rep(NA,n^2), nrow = nrow(h), ncol = ncol(h))
  for(i in seq_along(RHO)) {
    RHO[i] <- (1-alpha) * exp(-h[i]/a)
  }
  #diag(RHO) <- 1
  RHO
}
RHO <- get.RHO(MLIST.obs, h = h)

# SIGMA_xx (qN x qN)  q = number of soil properties = 9
SIGMA0.obs <- computeSigmaHat.LISREL(MLIST = MLIST.obs)[1:9,1:9] # qxq
SIGMA.xx <- kronecker(SIGMA0.obs, RHO)   # qN x qN
plotMat(SIGMA.xx[1:70,1:70])

# SIGMA_yy (qN x q)

SIGMA.yy <- computeSigmaHat.LISREL(MLIST = MLIST.obs)[1:9,1:9]

############## get prediction and residuals
# function to get prediction at new locations from trend model
# it includes parallel processing
get.pred <- function (MLIST = NULL, covar = covar.st){
  var.names <- c("CEC.Ar","CEC.Br","CEC.Cr","OC.Ar","OC.Br","OC.Cr",
                 "clay.Ar","clay.Br","clay.Cr","dem","vdchn","X",
                 "lstm","evisd","ndwi.b","twi","ndwi.a")
  m <- MLIST
  A <- m$beta[1:9,10:p]
  B <- m$beta[1:9,1:9]
  I <- diag(nrow = 9, ncol = 9)
  IB.inv <- solve(I - B)
  
  xy <- covar[,c("X","Y")]
  k <- covar[,var.names[10:p]]
  library(doParallel)
  registerDoParallel(6)
  pred <- foreach(i = icount(nrow(k)), .combine = rbind) %dopar%{
    p <- as.vector(as.matrix(k[i,]))
    t(IB.inv %*% A %*% p) 
  }
  doParallel::stopImplicitCluster()
  colnames(pred) <- var.names[1:9]
  cbind(pred,xy)
}
# function to get residuals
get.res <- function (m = NULL){
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

#var <- diag(SIGMA0.yy - SIGMA.yx %*% solve(SIGMA.xx) %*% SIGMA.xy)

# prediction with linear trend
pred <- get.pred(MLIST = MLIST.obs, covar = covar.st)
# residuals of the linear trend
res <- get.res(m = fit@Model@GLIST)

# SIGMA.xy (qN x q)
# at single location (latlong)
library(doParallel)
doParallel::registerDoParallel(cores = 3)
ks <- read.csv("Paper_4/data/ks.csv")[,c(colnames(s), "Y")] #
coord <- ks[,c("X", "Y")] 
y.all <- as.vector(res[,1:9]) # residuals SP

backup <- covar.st
loc <- covar.st[sample(x = rownames(covar.st),100), c("X","Y")]

# system.time({
k.res <-
  foreach(i = icount(nrow(loc)), .combine = rbind,
          .packages = c("sp")) %dopar%
 {
            ll <- loc[i,]
            coord.ll <- rbind(coord,ll)
            coord.ll$Y <- coord.ll$Y * ST$std.dev[21] / ST$std.dev[20]
            coordinates(coord.ll) <- ~X+Y
            # h = n x n matrix of distances between samples 
            h.all <- sp::spDists(coord.ll)
            h0 <- matrix(h.all[1:N,N+1], ncol = 1, nrow = N)
            #
            RHO0 <- get.RHO0(MLIST.obs, h = h0) # note that h0 is Nx1 
            SIGMA.xy <- kronecker(SIGMA.yy, RHO0) # qN x q
            SIGMA.yx <- t(SIGMA.xy) # q x qN
            t((SIGMA.yx %*% solve(SIGMA.xx) %*% y.all) )#* STt$std.dev[c(5:7,8:10,2:4)] )#+ 
            #STt$mean[c(5:7,8:10,2:4)])
          }
#})
colnames(k.res) <- colnames(res)
sd <- matrix(rep(STt$std.dev[c(5:7,8:10,2:4)], each=100), ncol=9)
means <- matrix(rep(STt$mean[c(5:7,8:10,2:4)], each=100), ncol=9)
summary((pred[rownames(loc),1:9] + k.res) *  sd + means)


(pred[rownames(loc),1:9] + predictions) * t(as.data.frame(STt$std.dev[c(5:7,8:10,2:4)]))
STt$mean[c(5:7,8:10,2:4)]


doParallel::stopImplicitCluster()

######################################################
####### Parameters bootstrapping #####################
######################################################
par <- read.csv("parameters.csv")[,-1]
#x <- rbind(as.data.frame(parameters),as.data.frame(x0))
library(lavaan)
partable <- partable(my.fit.lv.ML)
colnames(par) <- c(paste0(partable$lhs[partable$free!=0],
                          partable$op[partable$free!=0],
                          partable$rhs[partable$free!=0]), "alpha", "a")
# statistics of the estimates
library(reshape)
xmean <- as.data.frame(colMeans(par))
xmean <- melt(xmean)
xmean$variable <- names(par)

xsd <- data.frame(sd = apply(par, 2, sd))
xsd$min <-xmean$value - xsd$sd
xsd$max <-xmean$value + xsd$sd
xsd <- xsd[,-1]
xsd$variable <- rownames(xsd)

est <- as.data.frame(out$par)
est$variable <- names(par)
names(est) <- c("value", "variable")

lav <- data.frame(value = c(partable$est[partable$free!=0], NA, NA))
lav$variable <- names(par)



library(ggplot2)
meltx <- reshape::melt(par)

png(filename = "estimates.png",
    width = 2500, height = 2400, res = 180)
ggplot2::ggplot(data = meltx, mapping = aes(x = value)) +
  geom_rect(data = xsd, inherit.aes = F,
            aes(xmin=min, xmax = max, ymin = -Inf, ymax = +Inf), 
            fill = 'red', alpha = 0.2) +
  geom_histogram(bins = 40) + facet_wrap(~variable, scales = 'free_x') +
  geom_vline(aes(xintercept=value, color='Bootstrap'), data=xmean) + 
  geom_vline(aes(xintercept=value, color='Estimates'), data=est) +
  geom_vline(aes(xintercept=value, color='lavaan'), data=lav) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text=element_text(family="Serif"))
dev.off()

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
#   res <- matrix(data = NA, nrow = N, ncol = 9)
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


#### PERDICTION PROCESS ####
# comes from get.pred.points.w.covar.R 
# setwd("/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling/")
# 
# library(raster)
# #library(sp)
# 
# ## Points over DEM and its derivates
# files <- list.files(pattern = ".dat$")
# header <- gsub(".sdat", "", files)
# header <-  c("dem", "twi", "vdchn") 
# 
# #define crs
# wgs84 <- CRS("+init=epsg:4326")
# #UTM14N <- CRS("+init=epsg:32614")
# modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
# NAD83.KS.N <- CRS("+init=epsg:2796")
# 
# library(maptools)
# sa <- readShapeSpatial("/mnt/L0135974_DATA/UserData/BaseARG/study area/USA/Platte_area_extended.shp")
# 
# spplot(sa)
# proj4string(sa) <- wgs84
# sa <- spTransform(sa, NAD83.KS.N)
# extent(sa)
# 
# xcoord <- as.vector(extent(sa))[1:2]
# x <- seq(from=xcoord[1], to = xcoord[2], 1000)
# ycoord <- as.vector(extent(sa))[3:4]
# y <- seq(from = ycoord[1], to = ycoord[2], 1000)
# coord <- expand.grid(x,y)
# 
# xy <- SpatialPointsDataFrame(coords = coord,data =  as.data.frame(rep(NA, length(coord[,1]))))
# proj4string(xy) <- NAD83.KS.N
# 
# xy <- xy[rownames(over(sa, xy, returnList = TRUE)$`0`),]
# 
# 
# for(i in seq_along(files)){
#   xy@data[,i] <- extract(x = raster(files[i]), y = xy)
#   names(xy@data)[i] <- header[i]
# }
# 
# # list of tif files (modis)
# files <- list.files(pattern = ".tif$")[c(-1,-4)]
# header <- gsub(".tif", "", files)
# header <-  c("evisd", "lstm", "ndwi.a", "ndwi.b") 
# 
# # transform projection
# xy <- spTransform(x = xy, modis)
# 
# xy@data[,3+seq_along(files)] <- NA
# names(xy@data)[4:7] <- header
# 
# # extract values from tiff files
# for(i in seq_along(files)){
#   xy@data[,3+i] <- extract(x = raster(files[i]), y = xy)
# }
# xy <- spTransform(x = xy, NAD83.KS.N)
# 
# # r <- raster(xy, res = 1000)
# # proj4string(r)
# # plot(rasterize(xy[1],r, background= NA))
# # 
# # distanceFromPoints(object = , xy = , filename = )
# 
# xy <- as.data.frame(xy)
# names(xy)[8:9] <- c("X", "Y")
# xy$twi <- log10(xy$twi)
# xy$vdchn <- log10(xy$vdchn+10)
# xy$ndwi.a <- (xy$ndwi.a+10)^.3
# 
# #ST <- read.csv("~/Documents/SEM2DSM1/Paper_4/data/STt.ks.csv")
# 
# xy <- xy[,colnames(s)[10:18]]
# xy <- as.matrix(xy)
# #
# 
# #
# write.csv(file = "~/Documents/SEM2DSM1/Paper_4/data/xy.csv",x = xy)
