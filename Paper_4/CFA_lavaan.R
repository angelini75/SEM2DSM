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



# library(lavaan)
# HS.model <- ' visual  =~ x1 + x2 + x3
# textual =~ x4 + x5 + x6
# speed   =~ x7 + x8 + x9 
# textual ~ speed + visual
# speed ~ visual'
# 
# fit <- sem(HS.model, data=HolzingerSwineford1939, estimator = "ML") # or ML

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

x2MLIST(out.ML$par, MLIST)

# chi-square
out.ML$objective * 153
# 85.30551 
# library(DEoptim)
# results <- DEoptim(fn = objective_ML, lower = start.x-1, upper = start.x+1, 
#                    control = list(reltol=10E-15, steptol=50, itermax = 5000, trace = T, 
#                                   CR = 0.5, NP = 1000, F=0.8, strategy =2), MLIST = MLIST)
# chi-square
# results$optim$bestval * 153

summary(out.ML$par - results$optim$bestmem)

diff <- data.frame(par=1:51, diff=out.ML$par - results$optim$bestmem)


# > summary(out.ML$par - results$optim$bestmem)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.0175500 -0.0011810  0.0007089  0.0173700  0.0026310  0.4331000 

# par          diff
# par1    1  3.729837e-01
# par2    2 -1.754533e-02
# par3    3  4.331289e-01
# par4    4  8.171361e-02
# par5    5  5.895941e-04
# par6    6  1.160049e-03
# par7    7  1.957513e-03
# par8    8 -1.012679e-04
# par9    9 -5.287590e-03
# par10  10  2.166860e-03
# par11  11  1.149968e-03
# par12  12 -1.706141e-04
# par13  13  8.219489e-04
# par14  14 -3.769180e-04
# par15  15  4.117405e-03
# par16  16  1.118815e-03
# par17  17  2.083414e-03
# par18  18 -5.142453e-04
# par19  19 -1.294106e-03
# par20  20  1.916619e-03
# par21  21 -5.822167e-03
# par22  22  2.610429e-03
# par23  23  8.023352e-03
# par24  24 -2.596190e-04
# par25  25 -1.067525e-03
# par26  26  5.330471e-03
# par27  27  6.123758e-05
# par28  28  9.475898e-03
# par29  29  5.810359e-03
# par30  30  5.840557e-03
# par31  31  7.789799e-03
# par32  32 -7.485641e-03
# par33  33 -3.803355e-03
# par34  34  3.289150e-04
# par35  35  7.055161e-04
# par36  36 -5.983437e-03
# par37  37  7.088888e-04
# par38  38  2.651615e-03
# par39  39  2.004383e-03
# par40  40  7.173832e-04
# par41  41  4.783426e-03
# par42  42 -2.964310e-05
# par43  43  4.021224e-03
# par44  44 -7.375841e-03
# par45  45 -1.109567e-02
# par46  46  1.028914e-03
# par47  47  5.574870e-04
# par48  48 -1.659101e-03
# par49  49 -2.847471e-03
# par50  50 -3.047117e-03
# par51  51 -5.870095e-03