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
load("data_for_prediction_maps.RData")
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
  for(i in names(x)){
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
get.RHO <- function(MLIST = NULL, h = NULL) {
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
objective_ML <- function(x, MLIST = NULL) {
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
get.pred <- function (MLIST = NULL, covar = NULL){
  var.names <- c("CEC.Ar","CEC.Br","CEC.Cr","OC.Ar","OC.Br","OC.Cr",
                 "Clay.Ar","Clay.Br","Clay.Cr","dem","vdchn","X",
                 "lstm","evisd","ndwi.b","twi","ndwi.a")
  m <- MLIST
  A <- m$beta[1:9,10:p]
  B <- m$beta[1:9,1:9]
  I <- diag(nrow = 9, ncol = 9)
  IB.inv <- solve(I - B)
  k <- covar[,var.names[10:p]]
  pr <- as.vector(as.matrix(k))
  pred <- t(IB.inv %*% A %*% pr)
  colnames(pred) <- var.names[1:9]
  pred
}
# function to get residuals
get.res <- function (m = NULL, z = NULL){
  A <- m$beta[1:9,10:p]
  B <- m$beta[1:9,1:9]
  I <- diag(nrow = 9, ncol = 9)
  IB.inv <- solve(I - B)
  sp <- z[,1:9]
  pr <- z[,10:p]
  res <- matrix(data = NA, nrow = N, ncol = 9)
  for(i in seq_along(pr[,1])){
    res[i,] <- t(sp[i,] - (IB.inv %*% A %*% pr[i,]))
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

N = nrow(ks)
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
variance <- predicted
var.names <- c("CEC.Ar","CEC.Br","CEC.Cr","OC.Ar","OC.Br","OC.Cr",
               "Clay.Ar","Clay.Br","Clay.Cr","dem","vdchn","X",
               "lstm","evisd","ndwi.b","twi","ndwi.a")
h <- matrix()
z.all <- numeric()

cal <- ks[,]
rownames(cal) <- 1:N
# get h
xy <- cal[, c("X","Y")] # coordinates of cal data
xy[,"Y"] <- xy[,"Y"] * ST$std.dev[21] / ST$std.dev[20]
coordinates(xy) <- ~X+Y
h <- sp::spDists(xy) # h of cal data
# Optimize the function
cal <- as.matrix(cal[,colnames(s)])
z.all <- as.vector(cal)
out <- nlminb(start = start.x, objective = objective_ML,
              MLIST = MLIST, control = list(iter.max = 200)) 
MLIST.obs <- x2MLIST(out$par, MLIST) # extract parameters
RHO <- get.RHO(MLIST = MLIST.obs, h = h) # estimate RHO for the 146 samples
# # # Compute matrices: SIGMA.yy (9x9), SIGMA.xx (9Nx9N), SIGMA.xy (9x9N)
SIGMA.yy <- computeSigmaHat.LISREL(MLIST = MLIST.obs)[1:9,1:9] # qxq
SIGMA.xx <- kronecker(SIGMA.yy, RHO)   # qN x qN
#################################################################
################ prediction locations ##########################
###############################################################
s0.un <- read.csv("xy.csv")[,-1]
s0 <- std(x = s0.un, st = ST)
summary(s0)
xy.s0 <- s0[, c("X","Y")] # coordinates of cal data
xy.s0[,"Y"] <- xy.s0[,"Y"] * ST$std.dev[21] / ST$std.dev[20]


#################################################################
############### Foreach loop ###################################
###############################################################
library(doParallel)
registerDoParallel(cores = 12) # quics server, 48 cores: 6.5 hours

result <-
  foreach(i = icount(nrow(s0)), .combine = rbind, # 147 calibrations 
          .packages = "sp") %dopar% {
            # covar for prediction
            covar.st <- s0[i,c("dem","vdchn","X", "lstm","evisd","ndwi.b","twi","ndwi.a")]
            # prediction from linear model and residuals
            pred.lm <- get.pred(MLIST = MLIST.obs, covar = covar.st) 
            res <- get.res(m = MLIST.obs, z = s)
            y.all <- as.vector(res) # vector of residuals
            # Locations of calibration profiles (xy) and predition profile (ll)
            ll <- s0[i,c("X","Y")] # location of ks[i,]
            ll[,"Y"] <- ll[,"Y"] * ST$std.dev[21] / ST$std.dev[20]
            xy <- as.data.frame(xy)
            xy.ll <- rbind(xy,ll)
            coordinates(xy.ll) <- ~X+Y # all coordinates (N+1 = 147)
            h.all <- sp::spDists(xy.ll) 
            h0 <- matrix(h.all[1:N,N+1], ncol = 1, nrow = N) # h for prediction location
            RHO0 <- get.RHO0(MLIST.obs, h0 = h0) # note that h0 is Nx1
            # get SIGMA.xy
            SIGMA.xy <- kronecker(SIGMA.yy, RHO0) # qN x q
            # kriging of residuals
            k.res <- t(crossprod(SIGMA.xy, chol2inv(chol(SIGMA.xx))) %*% y.all) 
            colnames(k.res) <- paste0(colnames(res),".res")
            # total prediction
            predicted <- pred.lm + k.res
            predicted <- cbind(predicted, ll)
            # Computing prediction variance (could be in a function get.var())
            # PSI <- MLIST.obs$psi[1:9,1:9] # system error of SP
            # B <- MLIST.obs$beta[1:9,1:9] 
            # I <- diag(nrow = 9, ncol = 9)
            # IB.inv <- solve(I - B)
            # theta <- MLIST.obs$theta[1:9,1:9] # measurement error
            # var.SIGMA.yy <- PSI + IB.inv %*% tcrossprod(theta, IB.inv) #
            # var.SIGMA.xx <- kronecker(var.SIGMA.yy, RHO)
            # var.SIGMA.xy <- kronecker(var.SIGMA.yy, RHO0)
            # var.zeta <- var.SIGMA.yy -
            #   crossprod(var.SIGMA.xy,
            #             chol2inv(chol(var.SIGMA.xx))) %*% var.SIGMA.xy
            # variance <- matrix(diag(var.zeta), nrow = 1)
            # colnames(variance) <- paste0(colnames(res),".var")
            # #as.vector(var.zeta)
            result <- predicted #cbind(predicted, k.res, variance)
            result
          }

doParallel::stopImplicitCluster()
##########################################################################
################### Compute ME, RMSE and AVE ############################
########################################################################
result <- as.data.frame(result)
names(result)[1:9] <- gsub("r", "", names(result)[1:9])
write.csv(result, "prediction_spatSEM.csv")
result.sp <- result[,1:9]
xy <- result[,10:11]
ST <- ST[,-1]
# function to unstandardise the data
unstd<- function(x, st){
  y <- x
  for(i in names(x)){
    y[,i] <- (x[,i] * st[i,2]) + st[i,1]
  }
  y
}
# STsub <- ST[c(5:7,8:10,2:4),]

pred <- unstd(x = result.sp, st = ST)

xy.un <- (xy * ST["X",2])
xy.un$X <- xy.un$X + ST["X",1]
xy.un$Y <- xy.un$Y + ST["Y",1]

pred.xy <- cbind(pred,xy.un)

coordinates(pred.xy) <- ~X+Y
proj4string(pred.xy) <- CRS("+init=epsg:2796")

library(raster)
r <- raster(pred.xy, res = 1000.1)
proj4string(r) 
r.spatSEM <- rasterize(pred.xy,r, background= NA)

library(rasterVis)
levelplot(x = r.spatSEM[[c(5:7)]],
          ylab = "", xlab = "",
          names.attr= c("Clay of horizon A", "Clay of horizon B", "Clay of horizon C"),
          scales=list(draw=T),
          col.regions = colorRampPalette( c('#006d2c',"#edf8fb")),
          par.settings=list(grid.pars=list(fontfamily="serif")))

################################################################################
######################## Prediction SEM #######################################
##############################################################################

library(doParallel)
registerDoParallel(cores = 12) # quics server, 48 cores: 6.5 hours

result.SEM <-
  foreach(i = icount(nrow(s0)), .combine = rbind, # 147 calibrations 
          .packages = "sp") %dopar% {
            # covar for prediction
            covar.st <- s0[i,c("dem","vdchn","X", "lstm","evisd","ndwi.b","twi","ndwi.a")]
            # prediction from linear model and residuals
            pred.lm <- get.pred(MLIST = MLIST, covar = covar.st) 
            # total prediction
            pred.lm 
          }

doParallel::stopImplicitCluster()

result.SEM <- as.data.frame(result.SEM)
names(result.SEM)[1:9] <- gsub("r", "", names(result.SEM)[1:9])
write.csv(result.SEM, "prediction_SEM.csv")

pred.SEM <- unstd(x = result.SEM, st = ST)

pred.SEM.xy <- cbind(pred.SEM,xy.un)

coordinates(pred.SEM.xy) <- ~X+Y
proj4string(pred.SEM.xy) <- CRS("+init=epsg:2796")

r.SEM <- rasterize(pred.SEM.xy,r, background= NA)


################################################################################
############################ Plots MAPS #######################################
##############################################################################
CEC <- stack(r.SEM[[2]],r.spatSEM[[2]])
CEC.color <-  c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1',
                '#88419d','#810f7c','#4d004b')
CEC.plot <- levelplot(CEC, layout=c(1, 2), names.attr=c('Standard SEM', 'Spatial SEM'),
                      ylab = "", xlab = "", scales=list(draw=FALSE, alternating= FALSE),
                      at=seq(0, 40, length.out=40), 
                      par.strip.text=list(font=0.7),
                      par.settings=list(grid.pars=list(fontfamily="serif")),
                      col.regions = colorRampPalette(CEC.color), 
                      main = expression("CEC of A horizon"~~"/ cmol"[c]~~"kg"^{-1}))

OC <- stack(r.SEM[[5]],r.spatSEM[[5]])
OC.color <-  c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d',
               '#238443','#006837','#004529')
OC.plot <- levelplot(OC, layout=c(1, 2), names.attr=c('Standard SEM', 'Spatial SEM'), 
                     ylab = "", xlab = "",scales=list(draw=FALSE, alternating= FALSE),
                     at=seq(0, 4, length.out=40), 
                     par.strip.text=list(font=0.7),
                     par.settings=list(grid.pars=list(fontfamily="serif")),
                     col.regions = colorRampPalette(OC.color), 
                     main = expression("OC of A horizon / %"))

Clay <- stack(r.SEM[[8]],r.spatSEM[[8]])
Clay.color <-  c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a',
                '#e31a1c','#bd0026','#800026')
Clay.plot <- levelplot(Clay, layout=c(1, 2), names.attr=c('Standard SEM', 'Spatial SEM'), 
                       ylab = "", xlab = "", scales=list(draw=T, alternating= FALSE),
                       at=seq(0,50, length.out=50), 
                       par.strip.text=list(font=0.7),
                       par.settings=list(grid.pars=list(fontfamily="serif")),
                       col.regions = colorRampPalette(Clay.color), 
                       main = expression("Clay of A horizon / %"))

mat <- matrix(1:3, nrow=1)
plots <- list(CEC.plot,OC.plot,Clay.plot)

png(filename = "~/big/Fig5.png",
     width = 2500*1.40, height = 1270*1.40, res =  300)
for (i in 1:3){
  print(plots[[i]], split=c(col(mat)[i], row(mat)[i], ncol(mat), nrow(mat)), more=(i<3))
}
dev.off()

