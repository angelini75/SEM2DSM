############### #### ### ## # - SEM for fourth paper - # ## ### #### ###########
# Purpose        : To develope the code to calibrate a SE model with DEoptim
#                  taking into account spatial process on Zeta (system error)
# Maintainer     : Marcos Angelini  (marcos.angelini@wur.nl); 
# Contributions  : Gerard?
# Status         : 
# Note           : 
# sessionInfo(@RStudio desktop)  lenovo ThinkPad T430 (4 cores)
# R version 3.3.3 (2017-03-06)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 16.04.2 LTS
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
# [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] lavaan_0.5-23.1097
# 
# loaded via a namespace (and not attached):
#   [1] tools_3.3.3    mnormt_1.5-5   pbivnorm_0.6.0 stats4_3.3.3   quadprog_1.5-5

##### message from Yves at [google group]<https://groups.google.com/forum/#!searchin/lavaan/lavaan$20function$20github%7Csort:relevance/lavaan/0Hitqi3k_do/pYfyoocABwAJ>

#### How SEM models are estimated? ####

# The best way to learn it, is to program it yourself. Really. To get you
# started, I will paste below some code that uses some lavaan at the
# beginning (to set up the list of model matrices in MLIST), and to get
# some starting values, but apart from that, it is plain R. And we will
# fit the three-factor CFA model, using two different estimators: GLS and
# ML. Other estimators are easy to add. lavaan is nothing more but a
# slightly more fancy shell around this bare-bones code.
# 
# Hope this helps,
# 
# Yves.


# bare bones code to fit the 3-factor Holzinger & Swineford CFA model
# manually (YR 25 march 2016)
name <- function(x) { as.data.frame(names(x))}

library(lavaan)
HS.model <- ' visual  =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9 '

fit <- my.fit.lv.ML

fit <- cfa(HS.model, data=HolzingerSwineford1939, estimator = "ML") # or ML
partable(fit)
par.list <- inspect(fit)
# initial list of model matrices
MLIST <- lavTech(fit, "start")

# sample covariance matrix ####
S     <- fit@SampleStats@cov[[1]]

# inverse and logdet of S ####
# Choleski Decomposition
cS <- chol(S)
# inverse
S.inv <- chol2inv(cS)
# diagonal Choleski Decomposition
diag.cS <- diag(cS)
# logset 
S.logdet <- sum(log(diag.cS * diag.cS))

# Function to replace the sample -lambda, -theta and -psi by the estimated par. 
# x is the parameter vector: fill in new values in the model matrices
# par.list show the list of parameters
par.list
# free par. of beta
b.x <- par.list$beta[which(par.list$beta != 0)]
# free par. of psi
p.x <- as.vector(par.list$psi)[as.vector(par.list$psi)!=0]

x2MLIST <- function(x, MLIST) {
  beta.x <- x[b.x] # lambda: link between observed and latent variables
  #lambda.x <- x[0] # lambda: link between observed and latent variables
  #theta.x <- x[] # theta: variance of observed variables
  psi.x <- x[b.x] # var-covar matriz of latent var.
  
  #MLIST$lambda[ c(2,  3, 14, 15, 26, 27) ]            <- lambda.x
  #MLIST$theta[ c(1, 11, 21, 31, 41, 51, 61, 71, 81) ] <- theta.x
  MLIST$psi[which(as.vector(par.list$psi)!=0) ]           <- psi.x
  MLIST$beta[which(as.vector(par.list$beta)!=0)]          <- beta.x
  MLIST
}
# how to compute SigmaHat

get.IB.inv <- function (MLIST = NULL) {
  BETA <- MLIST$beta
  nr <- nrow(MLIST$psi)
  # if (!is.null(BETA)) {
  #   tmp <- -BETA
  #   tmp[lav_matrix_diag_idx(nr)] <- 1
  #   IB.inv <- solve(tmp)
  # }
  # else {
    IB.inv <- diag(nr)
  #}
  IB.inv
}


computeSigmaHat.LISREL <- function (MLIST = NULL, delta = TRUE) 
{
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  PSI <- MLIST$psi
  THETA <- MLIST$theta
  BETA <- MLIST$beta
  library(lavaan)
  if (is.null(BETA)) {
    LAMBDA..IB.inv <- LAMBDA
  }
  else {
    # lavaan:::.internal_get_IB.inv(MLIST = MLIST)
    IB.inv <- get.IB.inv(MLIST = MLIST)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }
  VYx <- tcrossprod(LAMBDA..IB.inv %*% PSI, LAMBDA..IB.inv) + 
    THETA
  if (delta && !is.null(MLIST$delta)) {
    DELTA <- diag(MLIST$delta[, 1L], nrow = nvar, ncol = nvar)
    VYx <- DELTA %*% VYx %*% DELTA
  }
  VYx
}
# objective function 'GLS'
# objective_GLS <- function(x, MLIST) {
#   
#   MLIST <- x2MLIST(x = x, MLIST = MLIST)
#   
#   # compute Sigma.hat
#   Sigma <- lavaan:::computeSigmaHat.LISREL(MLIST = MLIST)
#   
#   # GLS
#   R <- (S - Sigma)
#   A <- R %*% S.inv
#   objective <- 0.5 * sum(diag(A %*% A))
#   cat("objective = ", objective, "\n")
#   
#   objective
# }

# objective function 'ML'
objective_ML <- function(x, MLIST) {
  
  MLIST <- x2MLIST(x = x, MLIST = MLIST)
  
  # compute Sigma.hat
  Sigma <- lavaan:::computeSigmaHat.LISREL(MLIST = MLIST)
  nvar <- NROW(Sigma)
  
  # ML
  if (rcond(Sigma) < 0.003) {
# if (rcond(V) < 0.02) {# V is the matrix that is near singular, 
                        # you check the near singularity with the reciprocal 
                        # condition number and you define a threshold, in my case 0.02
    objective <- Inf  # objective is what you return at the end of you function,
                      # tell Inf and the parameters that lead to this near-singular 
                      # matrix are discarded.
    return(objective)
  }
  cS <- chol(Sigma)
  Sigma.inv <- chol2inv(cS)
  diag.cS <- diag(cS)
  logdet <- sum(log(diag.cS * diag.cS))
  objective <- logdet + sum(diag(S %*% Sigma.inv)) - S.logdet - nvar
  cat("objective = ", objective, "\n")
  
  objective
}

# get starting values
start.x <- parTable(fit)$start[ parTable(fit)$free > 0 ]

# # estimate parameters GLS
# out.GLS <- nlminb(start = start.x, objective = objective_GLS,
#                   MLIST = MLIST)
# out.GLS
# chi-square
#out.GLS$objective * 301
# 77.72896

# estimate parameters ML
out.ML  <- nlminb(start = start.x, objective = objective_ML,
                  MLIST = MLIST)

out.ML$par

library(DEoptim)
results <- DEoptim(fn = objective_ML, lower = start.x-1, upper = start.x+1, 
        control = list(reltol=10E-10, steptol=50, itermax = 5000, 
                       NP = 1000), MLIST = MLIST)
# trace = T, CR = 0.5, F = 0.8, NP = length(start.x)*10,
as.vector(results$optim$bestmem) - start.x

par1          par2          par3          par4          par5          par6          par7          par8 
0.8137945249  0.2275433953 -0.5394368627 -0.9855358009 -0.4572854421  0.1146149697  0.8188679047  0.6980586881 
par9         par10         par11         par12         par13         par14         par15         par16 
-0.1495758998 -0.2776271900 -0.2594496961 -0.0722540881  0.9862916521  0.0686602347  0.0362117604 -0.1123839995 
par17         par18         par19         par20         par21         par22         par23         par24 
0.1970414548  0.8486568404  0.7044257740  0.3136987149 -0.2651711782 -0.4004523634  0.0434787543  0.1397953289 
par25         par26         par27         par28         par29         par30         par31         par32 
-0.0001090689  0.5038536472  0.6166599369  0.0417478611 -0.8988876161  0.7344321572  0.2284543525 -0.7862457926 
par33         par34         par35         par36         par37         par38         par39         par40 
-0.1531693350  0.2468755269  0.5751980175  0.0762412744 -0.1219808293  0.6072973992  0.8519279457 -0.8720859538 
par41         par42         par43         par44         par45         par46         par47         par48 
-0.0317894349  0.6841117246 -0.8804179875 -0.1922467481 -0.7238852599  0.1114504316 -0.3029509712 -0.1824275297 
par49         par50         par51 
0.0290091560 -0.1577637694  0.2023927640 


out.ML$par
# chi-square
out.ML$objective * 301
# 85.30551 



###########
# Model 4 ####
# from Arg2Ks_6models.Rmd

#==============================================================================#
#### This code come from SEM_Arg2KS_1.0.R ####
#==============================================================================#

# Libraries ####
library(lavaan)
library(pastecs)
library(utils)
library(lattice)

rm(list=ls()[ls()!="t"])
name <- function(x) { as.data.frame(names(x))}
# chose one


# Dictionary of elements in this script ######
# d = calibration dataset. It comes from replacement_of_NAs.Rm 
#     (different versions: 5.0 to 5.3)
# ST = original mean and standard deviation of all variables
# STt = mean and standard deviation of transformed data
# nor = normalisation funcion (x-mean)/sd
# D = transformed and normalised data
# my.model = lavaan syntax
# my.fit = model fitted
# mod = modification indices (for respecification)
#------------------------------------------------#

e <- read.csv("~/Documents/SEM2DSM1/Paper_2/data/calib.data-5.0.csv")[,c(-1,-20)]
# Descriptive statistics and normality test. ####
round(stat.desc(e,norm = TRUE),3)
# Soil properties does not present strong deviation from normality.
# But some covariates need to be transformed. First, we store original mean and 
# sd in ST
ST.arg <- t(stat.desc(e,norm = TRUE)[c(9,13),])

# Based on normtest.W the following covariates need to be transformed
e$wdist <- e$wdist^0.5
e$maxc <- (e$maxc+20000)^2
e$slope <- e$slope^0.25
e$twi <- log10(e$twi)
e$vdchn <- log10(e$vdchn+10)
e$ndwi.a <- (e$ndwi.a+10)^.3
# New statistics
round(stat.desc(e,norm = TRUE),3)
# New mean and sd
STt.arg <- t(stat.desc(e,norm = TRUE)[c(9,13),])

# standardised data set ####
std <- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] - st[i,1]) / st[i,2]
  }
  y
}
Arg <- std(e,STt.arg)
Arg[,1] <- e[,1] 

setwd("~/Documents/SEM2DSM1/Paper_3/data/")
d <- read.csv("KS.data-0.2.csv")[,c(-1)] 
name(d)
names(d)[5:10] <- c("CEC.A","CEC.B","CEC.C","OC.A","OC.B","OC.C")
# remove outlayers
# d <- d[d$idp!=26058,]
# d <- d[d$idp!=22961,]
d <- cbind(d[1],
           d[,colnames(Arg)[2:10]],
           d[11:21])
d <- d[c(-156,-157),]
# Descriptive statistics and normality test. ####
round(stat.desc(d,norm = TRUE),3)
# Soil properties does not present strong deviation from normality.
# But some covariates need to be transformed. First, we store original mean and 
# sd in ST
ST.ks <- t(stat.desc(d,norm = TRUE)[c(9,13),])

# Based on normtest.W the following covariates need to be transformed
d$twi <- log10(d$twi)
d$vdchn <- log10(d$vdchn+10)
d$ndwi.a <- (d$ndwi.a+10)^.3
# OC as log10 of OC
# d$OC.A <- log10(d$OC.A)
# d$OC.B <- log10(d$OC.B)
# d$OC.C <- log10(d$OC.C)

# New statistics
round(stat.desc(d,norm = TRUE),3)
# New mean and sd
STt.ks <- t(stat.desc(d,norm = TRUE)[c(9,13),])

# standardised data set ####
std <- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] - st[i,1]) / st[i,2]
  }
  y
}
ks <- std(d,STt.ks)
ks[,1] <- d[,1] 
#######################################

# Model with latent variables ####
## Model (re)specification
my.model.lv <- '
# Measurement model (lamda and epsilon)
#--------------------#
CEC.Ar =~ 1*CEC.A
CEC.Br =~ 1*CEC.B
CEC.Cr =~ 1*CEC.C
OC.Ar =~ 1*OC.A
OC.Br =~ 1*OC.B
OC.Cr =~ 1*OC.C
clay.Ar =~ 1*clay.A
clay.Br =~ 1*clay.B
clay.Cr =~ 1*clay.C
## Measurement error #
# CEC.A ~~ 0.12 * CEC.A
# CEC.B ~~ 0.25 * CEC.B
# CEC.C ~~ 0.3 * CEC.C
OC.A ~~ 0.1 * OC.A
OC.B ~~ 0.1 * OC.B
OC.C ~~ 0.1 * OC.C
clay.A ~~ 0.05 * clay.A
clay.B ~~ 0.05 * clay.B
clay.C ~~ 0.05 *clay.C

#--------------------#
# Structural model (gamma and betta matrices)
#--------------------#
clay.Cr ~ dem + vdchn + X + lstm 
clay.Ar ~ clay.Cr + evisd + lstm + ndwi.b #+ Y 
clay.Br ~ clay.Ar + clay.Cr + vdchn + twi + ndwi.b + Y

OC.Ar ~ clay.Ar + evisd + lstm + ndwi.b 
OC.Br ~ OC.Ar + clay.Br + evisd + lstm + ndwi.a + vdchn
OC.Cr ~ OC.Br 

CEC.Ar ~ OC.Ar + clay.Ar 
CEC.Br ~ clay.Br + 0.01*OC.Br
CEC.Cr ~ clay.Cr + 0*OC.Cr

#------------------#
# Model error covariance (Psi)
#------------------#
CEC.Ar ~~ CEC.Br + CEC.Cr
CEC.Cr ~~ CEC.Br
#OC.Cr ~~ 0*CEC.Br + 0*CEC.Cr + 0*CEC.Ar 

#------------------#
# lavaan suggestions
#------------------#
clay.Br  ~    lstm
clay.Ar ~ twi

CEC.Ar ~~ clay.Br
OC.Ar ~   clay.Br
OC.Br ~   dem
OC.Cr ~~   clay.Br

CEC.Cr  ~  ndwi.a
CEC.Br  ~  ndwi.a + dem
CEC.Ar  ~    dem
#------------------#
'


# Model calibration ####
my.fit.lv.ML <- sem(model = my.model.lv,data = ks, meanstructure = FALSE, 
                    fixed.x = T)
inspect(my.fit.lv.ML,"cov.lv")
# Model evaluation ####
summary(my.fit.lv.ML, fit.measures=TRUE, rsquare = T)
mod.ks <- modindices(my.fit.lv.ML,sort. = T)
mod.ks[mod.ks$mi>7 & (mod.ks$op == "~"|mod.ks$op == "~~"),]# & mod.ks$lhs == "clay.Cr",] 


# Model validation ####
#### Prediction ####
pre <- cbind(ks[1,], matrix(nrow=1,ncol= 9, data = NA,
                            dimnames = list(NULL,paste0(names(ks)[2:10],".p"))))
a <- pre[-(1:nrow(pre)),c(2:10)] #observed
b <- pre[-(1:nrow(pre)),c(2:10)] #predicted
v <- pre[-(1:nrow(pre)),c(2:10)] #variance(s)
resids <- pre[-(1:nrow(pre)),c(2:10)] #residuals
theta <- pre[-(1:nrow(pre)),c(2:10)] # for mean and median
Var <- pre[-(1:nrow(pre)),c(2:10)] # model variance (constant)

# Loop: cal is calibration data, pre is prediction place
# Matrix dedinition (Section 3.3 2nd paper and Fig. 5) #

# Running Prediction @ i location #
# p is a matrix with the 10 external drivers
for(i in seq_along(ks[,1])){ 
  cal <- ks[,] # not LOO. Set [-i,] for LOO 
  pre[i,] <- ks[i,]
  # Fiting #
  my.fit.lv.ML <- sem(model = my.model.lv,data = cal, fixed.x = T,
                      estimator = "ML")
  
  # Matrix of Beta coefficients
  B <- inspect(my.fit.lv.ML, "est")$beta[1:9,1:9]
  # Identity matrix (Kappa coefficients)
  I <- diag(nrow = 9, ncol = 9)
  # Matrix of Gamma coefficients
  A <- inspect(my.fit.lv.ML, "est")$beta[1:9,10:18]
  # Matrix of Psi coefficients (model error variance-covariance)
  V <- inspect(my.fit.lv.ML, "est")$psi[1:9,1:9]
  # Matrix of measurement error (Epsylon)
  Th <- inspect(my.fit.lv.ML, "est")$theta[1:9,1:9]
  IB <- solve(I - B)
  
  # Running Prediction @ i location #
  # p is a matrix with the 10 external drivers
  p = as.vector(as.matrix(pre[i,colnames(A)])) # values of covariates ordered
  p = matrix(p, nrow = 9, ncol = 1)           # by lavaan sequence
  # prediction
  pre[i,22:30] = t(IB %*% A %*% p) # key equation
  # calculate standarised squared standard error
  ## theta is standarised squared standard error
  ## theta = ((observed-predicted)^2)/error variance=(standard error^2)
  a[i,] <- pre[i,c(2:10)] # observed values
  b[i,] <- pre[i,c(22:30)] # predicted values
  v <- diag(IB%*%V%*%t(IB)+Th) # error variance (it is not diagonal!)
  resids[i,] <- a[i,] - b[i,] # residuals
  theta[i,] <- (resids[i,]^2)/v # theta
  # Error variance #
  Var[i,] <- diag(IB %*% V %*% t(IB))
}

# Model variance

summary(Var)
Var <- apply(Var, MARGIN = 2, FUN = mean)

# function to unstandardise the data
unstd<- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] * st[i,2]) + st[i,1]
  }
  y
}

# Accuracy measures ####
# Residuals #
Res <- cbind(pre[,1], unstd(pre[,2:10], STt.ks[2:10,]), unstd(pre[,22:30],
                                                              STt.ks[2:10,]))
Res <-  cbind(pre[,1], pre[,2:10] - pre[,22:30])

(var.Res <- var(Res[2:10]))

(psi <- inspect(my.fit.lv.ML, "est")$psi[1:9,1:9])
var.Res - psi


Res
#################################################
### Autocorrelation in SEM residuals ####

names(Res)[1] <- "id.p"

R <- merge(x = Res, y = unique(d[,c(1,20,21)]), by.x = "id.p", by.y="idp")

library(sp)
library(gstat)
# R as geo
coordinates(R) <- ~X+Y

#define crs
wgs84 <- CRS("+init=epsg:4326")
UTM14N <- CRS("+init=epsg:32614")
modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
NAD83.KS.N <- CRS("+init=epsg:2796")

# Assign projection
proj4string(R) <- NAD83.KS.N

# check zero distance between profiles
zerodist(R, zero=0.0)
R <- remove.duplicates(R, zero = 0.0, remove.second = TRUE)

g <- list()
vg <- list()
vgm <- list()
par(mfrow = c(3, 3), pty = "s", mar=c(4,5,2,2), family="serif")
# CEC.A
g[[1]] <-  gstat(formula = CEC.A ~ 1, data = R)
vg[[1]] <- variogram(g[[1]], width = 7000, cutoff = 450000, cressie = TRUE)
plot(vg[[1]], plot.numbers = TRUE)

# # choose initial variogram model and plot:
vgm[[1]] <- vgm(nugget = 20,
                psill= 1,
                range=100000,
                model = "Exp")
vgm[[1]] <- fit.variogram(vg[[1]], vgm[[1]], fit.method = 6)
plot(vg[[1]], vgm[[1]], main = "CEC.A")
attr(vgm[[1]], "SSErr")

# CEC.B
g[[2]] <-  gstat(formula = CEC.B ~ 1, data = R)
vg[[2]] <- variogram(g[[2]], width = 7000, cutoff = 400000, cressie = TRUE)
plot(vg[[2]], plot.numbers = TRUE)

# # choose initial variogram model and plot:
vgm[[2]] <- vgm(nugget = 5,
                psill= 30,
                range=100000,
                model = "Exp")
vgm[[2]] <- fit.variogram(vg[[2]], vgm[[2]], fit.method = 7)
plot(vg[[2]], vgm[[2]], main = "CEC.B")
attr(vgm[[2]], "SSErr")

# CEC.C
g[[3]] <-  gstat(formula = CEC.C ~ 1, data = R)
vg[[3]] <- variogram(g[[3]], width = 7000, cutoff = 400000, cressie = TRUE)
plot(vg[[3]], plot.numbers = TRUE)

# # choose initial variogram model and plot:
vgm[[3]] <- vgm(nugget = 5,
                psill= 30,
                range=50000,
                model = "Sph")
vgm[[3]] <- fit.variogram(vg[[3]], vgm[[3]], fit.method = 7)
plot(vg[[3]], vgm[[3]], main = "CEC.C")
attr(vgm[[3]], "SSErr")
vgm[[3]]

# OC.A
g[[4]] <-  gstat(formula = OC.A ~ 1, data = R)
vg[[4]] <- variogram(g[[4]], width = 15000, cutoff = 400000, cressie = TRUE)
plot(vg[[4]], plot.numbers = TRUE)

vgm[[4]] <- vgm(nugget = 0.05,
                psill= 0.4,
                range=50000,
                model = "Exp")
vgm[[4]] <- fit.variogram(vg[[4]], vgm[[4]], fit.method = 7)
plot(vg[[4]], vgm[[4]], main = "OC.A")
attr(vgm[[4]], "SSErr")
vgm[[4]]

# OC.B
g[[5]] <-  gstat(formula = OC.B ~ 1, data = R)
vg[[5]] <- variogram(g[[5]], width = 15000, cutoff = 400000, cressie = TRUE)
plot(vg[[5]], plot.numbers = TRUE)

vgm[[5]] <- vgm(nugget = 0,
                psill= 0.8,
                range=20000,
                model = "Exp")
vgm[[5]] <- fit.variogram(vg[[5]], vgm[[5]], fit.method = 7)
plot(vg[[5]], vgm[[5]], main = "OC.B")
attr(vgm[[5]], "SSErr")
vgm[[5]]


# OC.C
g[[6]] <-  gstat(formula = OC.C ~ 1, data = R)
vg[[6]] <- variogram(g[[6]], width = 10000, cutoff = 400000, cressie = TRUE)
plot(vg[[6]], plot.numbers = TRUE)

vgm[[6]] <- vgm(nugget = 0.1,
                psill= 0.5,
                range=20000,
                model = "Gau")
vgm[[6]] <- fit.variogram(vg[[6]], vgm[[6]], fit.method = 2)
plot(vg[[6]], vgm[[6]], main = "OC.C")
attr(vgm[[6]], "SSErr")
vgm[[6]]


# Clay.A
g[[7]] <-  gstat(formula = clay.A ~ 1, data = R)
vg[[7]] <- variogram(g[[7]], width = 7000, cutoff = 400000, cressie = TRUE)
plot(vg[[7]], plot.numbers = TRUE)

vgm[[7]] <- vgm(nugget = 20,
                psill= 80,
                range=100000,
                model = "Exp")
vgm[[7]] <- fit.variogram(vg[[7]], vgm[[7]], fit.method = 1)
plot(vg[[7]], vgm[[7]], main = "Clay.A")
attr(vgm[[7]], "SSErr")
vgm[[7]]

# Clay.B
g[[8]] <-  gstat(formula = clay.B ~ 1, data = R)
vg[[8]] <- variogram(g[[8]], width = 7000, cutoff = 400000, cressie = TRUE)
plot(vg[[8]], plot.numbers = TRUE)

vgm[[8]] <- vgm(nugget = 15,
                psill= 120,
                range=100000,
                model = "Sph")
vgm[[8]] <- fit.variogram(vg[[8]], vgm[[8]], fit.method = 1)
plot(vg[[8]], vgm[[8]], main = "Clay.B")
attr(vgm[[8]], "SSErr")
vgm[[8]]

# Clay.C
g[[9]] <-  gstat(formula = clay.C ~ 1, data = R)
vg[[9]] <- variogram(g[[9]], width = 15000, cutoff = 400000, cressie = TRUE)
plot(vg[[9]], plot.numbers = TRUE)

vgm[[9]] <- vgm(nugget = 20,
                psill= 100,
                range=100000,
                model = "Gau")
vgm[[9]] <- fit.variogram(vg[[9]], vgm[[9]], fit.method = 2)
plot(vg[[9]], vgm[[9]], main = "Clay.C")
attr(vgm[[9]], "SSErr")
vgm[[9]]


# Three graphs (soil properties)
# tiff(filename = "~/Dropbox/PhD Marcos/Paper 4/Figures/Fig1.tif", 
#       width = 2500, height = 1000, res =  350)
par(mfrow = c(1, 3), pty = "s", mar=c(4,5,2,2), family="serif")
### CEC
## A
plot(variogramLine(vgm[[1]], maxdist=500000), 
     type="l", lwd=2,col="#AA0000",
     main= "CEC", 
     xlab = "Distance / m", 
     ylab = expression("Semivariance"~~"/"~("cmol"[c]~~"kg"^{-1})^{2}),
     cex.lab = 1.3, ylim=c(0,1.2))
points(gamma ~ dist, vg[[1]], col="#770000")
#legend(x= "topleft",legend = "SSErr", bty = "n")
#legend(x= 12,legend = round(attr(vgm[[1]], "SSErr"),3), text.col ="#AA0000", 
#        bty = "n")
## B
lines(variogramLine(vgm[[2]], maxdist=500000), lwd=2, col="#00AA00")
points(gamma ~ dist, vg[[2]], col="#007700")
## C
lines(variogramLine(vgm[[3]], maxdist=500000), lwd=2, col="#0000AA")
points(gamma ~ dist, vg[[3]], col="#000077")

### OC
## A
plot(variogramLine(vgm[[4]], maxdist=500000), type="l", lwd=2,col="#AA0000", 
     main= "OC",
     xlab = "Distance / m", 
     ylab = expression("Semivariance"~~"/ %"^{2}),
     cex.lab = 1.3, ylim=c(0,1.2)) 
points(gamma ~ dist, vg[[4]], col="#770000")
## B
lines(variogramLine(vgm[[5]], maxdist=500000), lwd=2, col="#00AA00")
points(gamma ~ dist, vg[[5]], col="#007700")
## C
lines(variogramLine(vgm[[6]], maxdist=500000), lwd=2, col="#0000AA")
points(gamma ~ dist, vg[[6]], col="#000077")

### Clay
## A
plot(variogramLine(vgm[[7]], maxdist=500000), type="l", lwd=2,col="#AA0000", 
     main="Clay",
     xlab = "Distance / m", 
     ylab = expression("Semivariance"~~"/ %"^{2}),
     cex.lab = 1.3, ylim=c(0,1.2)) 
points(gamma ~ dist, vg[[7]], col="#770000")
## B
lines(variogramLine(vgm[[8]], maxdist=500000), lwd=2, col="#00AA00")
points(gamma ~ dist, vg[[8]], col="#007700")
## C
lines(variogramLine(vgm[[9]], maxdist=500000), lwd=2, col="#0000AA")
points(gamma ~ dist, vg[[9]], col="#000077")

#dev.off()

###### Cross-variograms ######

g <- list()
vg <- list()
vgm <- list()

# CEC.A
g[[10]] <-  gstat(id = c("CEC.A", "CEC.B", "CEC.C"), data = R)
vg[[1]] <- variogram(g[[1]], width = 7000, cutoff = 450000, cressie = TRUE)
plot(vg[[1]], plot.numbers = TRUE)

# # choose initial variogram model and plot:
vgm[[1]] <- vgm(nugget = 20,
                psill= 1,
                range=100000,
                model = "Exp")
vgm[[1]] <- fit.variogram(vg[[1]], vgm[[1]], fit.method = 6)
plot(vg[[1]], vgm[[1]], main = "CEC.A")

names(R)[8:10] <- c("Clay.A", "Clay.B", "Clay.C") 
rm(cv)
cv <- gstat(id = "CEC.A", formula = CEC.A ~ + 1, data = R, nmax = 10)
cv <- gstat(cv, id = "CEC.B", formula = CEC.B ~ + 1, data = R, nmax = 10)
cv <- gstat(cv, id = "CEC.C", formula = CEC.C ~ + 1, data = R, nmax = 10)
cv <- gstat(cv, id = "OC.A", formula = OC.A ~ + 1, data = R, nmax = 10)
cv <- gstat(cv, id = "OC.B", formula = OC.B ~ + 1, data = R, nmax = 10)
cv <- gstat(cv, id = "OC.C", formula = OC.C ~ + 1, data = R, nmax = 10)
cv <- gstat(cv, id = "Clay.A", formula = Clay.A ~ + 1, data = R, nmax = 10)
cv <- gstat(cv, id = "Clay.B", formula = Clay.B ~ + 1, data = R, nmax = 10)
cv <- gstat(cv, id = "Clay.C", formula = Clay.C ~ 1, data = R, nmax = 10)
cv <- gstat(cv, 
            model = vgm(nugget = 20,
                        psill= 1,
                        range=100000,
                        model = "Exp"), 
            fill.all = T)
cv
cv.var<- variogram(cv, cutoff = 450000) 
plot(cv.var)
names(meuse.g)
cv.fit<-fit.lmc(cv.var, cv) 

png(filename = "~/Dropbox/PhD Marcos/Paper 4/Figures/Fig2.png", 
     width = 3000, height = 3000, res =  250)
plot(cv.var, model=cv.fit, 
     main="Variograms and cross-variograms of standardized residuals",
     xlab = "Distance / m", 
     ylab = "Semivariance",
     scales=list(x = list(alternating = 1), y = list(alternating = 1)),
     par.settings=list(grid.pars=list(fontfamily="serif")))
dev.off()

