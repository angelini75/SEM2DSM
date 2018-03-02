# Get P-values of estimates for paper 4
rm(list=ls())
name <- function(x) { as.data.frame(names(x))}

# from http://www.cyclismo.org/tutorial/R/pValues.html ##############

# table with coefficients. est = standard sem; new.est = spatial sem
coef <- read.csv("~/git/SEM2DSM/Paper_4/data/est_se.csv")[,-1]
name(coef)
coef <- coef[c(1:6,8)]
# sem = standard sem; spt = spatial sem; se = standard error
names(coef) <- c("lhs","op","rhs","est.sem","se.sem","est.spt","se.spt")
coef$p.sem <- 0
coef$p.spt <- 0
coef <- coef[c(1:5,8,6,7,9)]

# P-value for sem
a <- coef$est.sem
s <- coef$se.sem
n <- 147
xbar <- 0
z <- (xbar-a)/s
coef$p.sem <- 2*pnorm(-abs(z))

# P-value for spt
a <- coef$est.spt
s <- coef$se.spt
n <- 147
xbar <- 0
z <- (xbar-a)/s
coef$p.spt <- 2*pnorm(-abs(z))



write.csv(coef, "~/git/SEM2DSM/Paper_4/data/est_se_pvalue.csv")
