# Get P-values of estimates for paper 4
# from http://www.cyclismo.org/tutorial/R/pValues.html

# table with coefficients. est = standard sem; new.est = spatial sem
coef <- read.csv("~/git/SEM2DSM/Paper_4/data/est_se.csv")
a <- coef$new.est
s <- coef$lav.se
n <- 147
xbar <- 0
z <- (xbar-a)/(s/sqrt(n))
z

p_value <- 2*pnorm(-abs(z))

coef$P_value <- p_value
coef$signif <- coef$P_value<0.05

write.csv(coef, "~/git/SEM2DSM/Paper_4/data/est_se_pvalue.csv")
