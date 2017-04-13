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



library(lavaan)
HS.model <- ' visual  =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9 
textual ~ speed + visual
speed ~ visual'

fit <- sem(HS.model, data=HolzingerSwineford1939, estimator = "ML") # or ML

fit <- my.fit.lv.ML
# initial list of model matrices
MLIST <- lavTech(fit, "start")

# sample covariance matrix
S     <- fit@SampleStats@cov[[1]]

# inverse and logdet of S
cS <- chol(S)
S.inv <- chol2inv(cS)
diag.cS <- diag(cS)
S.logdet <- sum(log(diag.cS * diag.cS))

# x is the parameter vector: fill in new values in the model matrices
par.list <- inspect(fit)
x2MLIST <- function(x, MLIST) {
  lambda.x <- x[as.vector(par.list$lambda)[as.vector(par.list$lambda)!=0]]
  theta.x <- x[as.vector(par.list$theta)[as.vector(par.list$theta)!=0]]
  psi.x <- x[as.vector(par.list$psi)[as.vector(par.list$psi)!=0]]
  beta.x <- x[as.vector(par.list$beta)[as.vector(par.list$beta)!=0]]
  
  MLIST$lambda[which(as.vector(par.list$lambda)!=0)]     <- lambda.x
  MLIST$theta[which(as.vector(par.list$theta)!=0)]       <- theta.x
  MLIST$psi[which(as.vector(par.list$psi)!=0)]           <- psi.x
  MLIST$beta[which(as.vector(par.list$beta)!=0)]           <- beta.x
  MLIST
}
# function
get.IB.inv <- function (MLIST = NULL) {
  BETA <- MLIST$beta
  nr <- nrow(MLIST$psi)
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
    lavaan:::.internal_get_IB.inv(MLIST = MLIST)
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

# objective function 'ML'

# Sigma <- matrix(data = c(0.59156973, 0.46014155, 0.65501854, 0.05000000, 0.05664468, 0.04620913, 0.05000000, 0.06125368, 0.04272028, 0.46014155, 0.35791257, 0.50949402, 0.03889157, 0.04406002, 0.03594292, 0.03889157, 0.04764504, 0.03322918, 0.65501854, 0.50949402, 0.72527254, 0.05536275, 0.06272011, 0.05116529, 0.05536275, 0.06782344, 0.04730225, 0.05000000, 0.03889157, 0.05536275, 1.19017633, 0.57889891, 0.47224938, 0.05000000, 0.06125368, 0.04272028, 0.05664468, 0.04406002, 0.06272011, 0.57889891, 1.34672282, 0.53500831, 0.05664468, 0.06939390, 0.04839754, 0.04620913, 0.03594292, 0.05116529, 0.47224938, 0.53500831, 1.07387710, 0.04620913, 0.05660959, 0.03948134, 0.05000000, 0.03889157, 0.05536275, 0.05000000, 0.05664468, 0.04620913, 1.18283419, 0.62172721, 0.43361254, 0.06125368, 0.04764504, 0.06782344, 0.06125368, 0.06939390, 0.05660959, 0.62172721, 1.59155446, 0.53120726, 0.04272028, 0.03322918, 0.04730225, 0.04272028, 0.04839754, 0.03948134, 0.43361254, 0.53120726, 0.96866021),
#                 nrow = 9, ncol = 9)
objective_ML <- function(x, MLIST) {
  MLIST <- x2MLIST(x = x, MLIST = MLIST)
  # compute Sigma.hat
  Sigma <- computeSigmaHat.LISREL(MLIST = MLIST)
  if (all(eigen(Sigma)$values >0)) {
    nvar <- NROW(Sigma)
    cS <- chol(Sigma)
    Sigma.inv <- chol2inv(cS)
    diag.cS <- diag(cS)
    logdet <- sum(log(diag.cS * diag.cS))
    objective <- logdet + sum(diag(S %*% Sigma.inv)) - S.logdet - nvar
    cat("objective = ", objective, "\n")
    objective
  } else {
    objective <- Inf
    objective
  }
}
# get starting values
start.x <- parTable(fit)$start[ parTable(fit)$free > 0 ]

# estimate parameters ML
out.ML  <- nlminb(start = start.x, objective = objective_ML, 
                  MLIST = MLIST, control = list(eval.max = 5000, trace = 1))

# chi-square
out.ML$objective * 156
# 85.30551 
library(DEoptim)
results <- DEoptim(fn = objective_ML, lower = start.x-1, upper = start.x+1, 
                   control = list(reltol=10E-15, steptol=50, itermax = 5000, trace = T, 
                                  CR = 0.5, NP = 1000, F=0.8, strategy =2), MLIST = MLIST)
# chi-square
results$optim$bestval * 156

summary(out.ML$par - results$optim$bestmem)

# > summary(out.ML$par - results$optim$bestmem)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.0175500 -0.0011810  0.0007089  0.0173700  0.0026310  0.4331000 