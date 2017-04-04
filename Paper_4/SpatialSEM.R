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
x2MLIST <- function(x, MLIST) {
  lambda.x <- x[1:6] # lambda: link between observed and latent variables
  theta.x <- x[7:15] # theta: variance of observed variables
  psi.x <- x[ c(16,19,20,19,17,21,20,21, 18) ] # var-covar matriz of latent var.
  
  MLIST$lambda[ c(2,  3, 14, 15, 26, 27) ]            <- lambda.x
  MLIST$theta[ c(1, 11, 21, 31, 41, 51, 61, 71, 81) ] <- theta.x
  MLIST$psi[ c(1, 2, 3, 4, 5, 6, 7, 8, 9) ]           <- psi.x
  
  MLIST
}
# how to compute SigmaHat
computeSigmaHat.LISREL <- function (MLIST = NULL, delta = TRUE) 
{
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  PSI <- MLIST$psi
  THETA <- MLIST$theta
  BETA <- MLIST$beta
  if (is.null(BETA)) {
    LAMBDA..IB.inv <- LAMBDA
  }
  else {
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
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
objective_GLS <- function(x, MLIST) {
  
  MLIST <- x2MLIST(x = x, MLIST = MLIST)
  
  # compute Sigma.hat
  Sigma <- lavaan:::computeSigmaHat.LISREL(MLIST = MLIST)
  
  # GLS
  R <- (S - Sigma)
  A <- R %*% S.inv
  objective <- 0.5 * sum(diag(A %*% A))
  cat("objective = ", objective, "\n")
  
  objective
}

# objective function 'ML'
objective_ML <- function(x, MLIST) {
  
  MLIST <- x2MLIST(x = x, MLIST = MLIST)
  
  # compute Sigma.hat
  Sigma <- lavaan:::computeSigmaHat.LISREL(MLIST = MLIST)
  nvar <- NROW(Sigma)
  
  # ML
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

# estimate parameters GLS
out.GLS <- nlminb(start = start.x, objective = objective_GLS,
                  MLIST = MLIST)
out.GLS
# chi-square
out.GLS$objective * 301
# 77.72896

# estimate parameters ML
out.ML  <- nlminb(start = start.x, objective = objective_ML,
                  MLIST = MLIST)
out.ML
# chi-square
out.ML$objective * 301
# 85.30551 