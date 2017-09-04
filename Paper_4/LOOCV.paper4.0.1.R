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
variance <- predicted
var.names <- c("CEC.Ar","CEC.Br","CEC.Cr","OC.Ar","OC.Br","OC.Cr",
               "Clay.Ar","Clay.Br","Clay.Cr","dem","vdchn","X",
               "lstm","evisd","ndwi.b","twi","ndwi.a")
h <- matrix()
z.all <- numeric()
#################################################################
############### Foreach loop ###################################
###############################################################
library(doParallel)
registerDoParallel(cores = 12) # quics server, 48 cores: 6.5 hours

result <-
  foreach(i = icount(nrow(ks)), .combine = rbind, # 147 calibrations 
          .packages = "sp") %dopar% {
            # calibration dataset
            cal <- ks[-i,]
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
                          MLIST = MLIST, control = list(iter.max = 200)) # ~150 iterations
            MLIST.obs <- x2MLIST(out$par, MLIST) # extract parameters
            RHO <- get.RHO(MLIST = MLIST.obs, h = h) # estimate RHO for the 146 samples
            # # # Compute matrices: SIGMA.yy (9x9), SIGMA.xx (9Nx9N), SIGMA.xy (9x9N)
            SIGMA.yy <- computeSigmaHat.LISREL(MLIST = MLIST.obs)[1:9,1:9] # qxq
            SIGMA.xx <- kronecker(SIGMA.yy, RHO)   # qN x qN
            # covar for prediction
            covar.st <- ks[i,c("dem","vdchn","X", "lstm","evisd","ndwi.b","twi","ndwi.a")]
            # prediction from linear model and residuals
            pred.lm <- get.pred(MLIST = MLIST.obs, covar = covar.st) 
            res <- get.res(m = MLIST.obs, z = cal)
            y.all <- as.vector(res) # vector of residuals
            # Locations of calibration profiles (xy) and predition profile (ll)
            ll <- ks[i,c("X","Y")] # location of ks[i,]
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
            # Computing prediction variance (could be in a function get.var())
            PSI <- MLIST.obs$psi[1:9,1:9] # system error of SP
            B <- MLIST.obs$beta[1:9,1:9] 
            I <- diag(nrow = 9, ncol = 9)
            IB.inv <- solve(I - B)
            theta <- MLIST.obs$theta[1:9,1:9] # measurement error
            var.SIGMA.yy <- PSI + IB.inv %*% tcrossprod(theta, IB.inv) #
            var.SIGMA.xx <- kronecker(var.SIGMA.yy, RHO)
            var.SIGMA.xy <- kronecker(var.SIGMA.yy, RHO0)
            var.zeta <- var.SIGMA.yy -
              crossprod(var.SIGMA.xy,
                        chol2inv(chol(var.SIGMA.xx))) %*% var.SIGMA.xy
            variance <- matrix(diag(var.zeta), nrow = 1)
            colnames(variance) <- paste0(colnames(res),".var")
            #as.vector(var.zeta)
            result <- cbind(predicted, k.res, variance)
            result
          }

doParallel::stopImplicitCluster()
##########################################################################
################### Compute ME, RMSE and AVE ############################
########################################################################
result <- as.data.frame(result)
write.csv(result, "results_LOOCV_2.csv")
#result <- read.csv("Paper_4/data/results_LOOCV.csv")[,-1]
names(result)[1:9] <- gsub("r", "", names(result)[1:9])
result.sp <- result[,1:9]#cbind(,ks[,c("X","Y")])
ST <- ST[,-1]
# function to unstandardise the data
unstd<- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] * st[i,2]) + st[i,1]
  }
  y
}
STsub <- ST[c(5:7,8:10,2:4),]

obs <- unstd(x = ks[,names(result.sp)], STsub)
pred <- unstd(x = result.sp, STsub)


# create report
report <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, SS = NA,
                     mean_theta = NA, median_th = NA)
obs <- unstd(x = ks[,names(result.sp)], STsub)
pred <- unstd(x = result.sp, STsub)
Res <- cbind(obs, pred)
names(Res)[10:18] <- paste0(names(Res)[10:18],".p")
for (i in 1:9) {
  # ME <- mean error 
  ME  <-  mean(Res[,i] - Res[,i + 9])
  # RMSE (root mean squared error)
  RMSE <- sqrt(mean((Res[,i] - Res[,i + 9]) ^ 2))
  MSE <- mean((Res[,i] - Res[,i + 9]) ^ 2)
  # SS (Sum of squares)
  SS <- sum((Res[,i] - Res[,i + 9]) ^ 2)
  # fill report table
  report[i,1] <- c(names(Res)[i])
  report[i,2] <- ME
  report[i,3] <- RMSE
  report[i,4] <- SS
}

# for(i in 1:9){
#   report$mean_theta[i] <- mean(theta[,i])
#   report$median_th[i] <- median(theta[,i])
# }
report
#d.stat <- read.csv("summary.calibdata.csv")
STsub$SS <- NA 
for(i in 1:9){
  STsub$SS[i] <- sum((obs[,i] - STsub$mean[i])^2)
}

report$R2 <- 1 - (as.numeric(report$SS) / as.numeric(STsub$SS))
report

# Analysis by Soil Property
# plot mesured vs predicted combined ####par(mfrow = c(1,3), pty="s",mai=rep(0.7,4))
par(mfrow = c(1, 3), pty="s",mai=rep(0.7,4))
rsq<- NULL
CEC <- rbind(as.matrix(Res[,c(1,10)]), as.matrix(Res[,c(2,11)]),
             as.matrix(Res[,c(3,12)]))
colnames(CEC) <- c("CECo","CECp")
rownames(CEC) <- 1:nrow(CEC)
CEC <- as.data.frame(CEC)
rsq[1] <- 1 - (sum((CEC$CECo - CEC$CECp)^2)/
                 sum((mean(CEC$CECo)-CEC$CECo)^2))

OC <- rbind(as.matrix(Res[,c(4,13)]), as.matrix(Res[,c(5,14)]),
            as.matrix(Res[,c(6,15)]))
colnames(OC) <- c("OCo","OCp")
rownames(OC) <- 1:nrow(OC)
OC <- as.data.frame(OC)
rsq[2] <- 1 - (sum((OC$OCo - OC$OCp)^2)/
                 sum((mean(OC$OCo)-OC$OCo)^2))

clay <- rbind(as.matrix(Res[,c(7,16)]), as.matrix(Res[,c(8,17)]),
              as.matrix(Res[,c(9,18)]))
colnames(clay) <- c("clayo","clayp")
rownames(clay) <- 1:nrow(clay)
clay <- as.data.frame(clay)
rsq[3] <- 1 - (sum((clay$clayo - clay$clayp)^2)/
                 sum((mean(clay$clayo)-clay$clayo)^2))


# create report by soil property
report2 <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, r2 = NA)
z <- cbind(CEC, OC, clay)
for (i in c(1,3,5)) {
  # ME <- mean error 
  ME  <-  mean(z[,i] - z[,i + 1])
  # RMSE (root mean squared error)
  RMSE <- sqrt(mean((z[,i] - z[,i + 1]) ^ 2))
  MSE <- mean((z[,i] - z[,i + 1]) ^ 2)
  # fill report table
  report2[i,1] <- c(names(z)[i])
  report2[i,2:3] <- c(ME, RMSE)
}

report2$r2[1] <- rsq[1]
report2$r2[3] <- rsq[2]
report2$r2[5] <- rsq[3]

report2 <- report2[c(-4,-2),]
names(report)[7] <- "r2"
report.total <- rbind(report[,c(1:3,7)],report2)

write.csv(report.total, "Paper_4/data/report_LOOCV_2.csv")

# kriged residuals statistics

kres <- result[, 10:18]

kres.un <- unstd.b(x = kres,st = STsub)
unstd.b<- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] * st[i,2])
  }
  y
}
round(colMeans(kres.un), 3)
round(apply(kres.un, FUN = sd, MARGIN = 2), 3)

################################################################################
# Plot with lattice ####
head(Res)
Res$X <- NA
Res <- Res[,c(19,1:18)]
residuales <- rbind(data.frame(sp="CEC", hor="Joint h.", Obs=Res[,2], Pred=Res[,11]),
                    data.frame(sp="CEC", hor="Joint h.", Obs=Res[,3], Pred=Res[,12]),
                    data.frame(sp="CEC", hor="Joint h.", Obs=Res[,4], Pred=Res[,13]),
                    data.frame(sp="OC", hor="Joint h.", Obs=Res[,5], Pred=Res[,14]),
                    data.frame(sp="OC", hor="Joint h.", Obs=Res[,6], Pred=Res[,15]),
                    data.frame(sp="OC", hor="Joint h.", Obs=Res[,7], Pred=Res[,16]),
                    data.frame(sp="Clay", hor="Joint h.", Obs=Res[,8], Pred=Res[,17]),
                    data.frame(sp="Clay", hor="Joint h.", Obs=Res[,9], Pred=Res[,18]),
                    data.frame(sp="Clay", hor="Joint h.", Obs=Res[,10], Pred=Res[,19]),
                    data.frame(sp="CEC", hor="C", Obs=Res[,4], Pred=Res[,13]),
                    data.frame(sp="CEC", hor="B", Obs=Res[,3], Pred=Res[,12]),
                    data.frame(sp="CEC", hor="A", Obs=Res[,2], Pred=Res[,11]),
                    data.frame(sp="OC", hor="C", Obs=Res[,7], Pred=Res[,16]),
                    data.frame(sp="OC", hor="B", Obs=Res[,6], Pred=Res[,15]),
                    data.frame(sp="OC", hor="A", Obs=Res[,5], Pred=Res[,14]),
                    data.frame(sp="Clay", hor="C", Obs=Res[,10], Pred=Res[,19]),
                    data.frame(sp="Clay", hor="B", Obs=Res[,9], Pred=Res[,18]),
                    data.frame(sp="Clay", hor="A", Obs=Res[,8], Pred=Res[,17]))

library(lattice)
library(latticeExtra)
library(hexbin)
png(filename = "~/big/Scatterplot_res.png", 
     width = 2000, height = 2000, res =  300)
plot <- xyplot(Pred ~ Obs| sp + hor, data=residuales, 
               type=c('p', "g"),asp = 0.9,.aspect.ratio = 1, 
               default.scales = list(tick.number=3, tck = 1, minlength = 3),
               scales = list(alternating= FALSE,
                             y=list(relation='free'), x=list(relation='free'),
                             limits=rep(list(c(0,40), c(0,4), c(0,70)), 4)),
               par.settings=list(grid.pars=list(fontfamily="serif")),
               pch = ".", cex = 3, alpha = 0.4, col = "black",
               xlab = "Observed", ylab = "Predicted",
               
)
#panel = panel.hexbinplot())

useOuterStrips(combineLimits(x = plot,
                             margin.x = c(2), margin.y = c(),
                             extend = T, adjust.labels = T),
               strip = strip.custom(
                 factor.levels = c(expression("CEC"~~"/ cmol"[c]~~"kg"^{-1}),
                                   "OC / %",
                                   "Clay / %")))
dev.off()
