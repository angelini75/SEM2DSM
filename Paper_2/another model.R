
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration")
rm(list=ls())
#install.packages('reshape')
#install.packages("lavaan")#, repos="http://www.da.ugent.be", type="source")
# install.packages('soiltexture')
#install.packages('corrgram')
# install.packages("plyr")
#library(soiltexture)
# library(reshape)
# library(plyr)
# library(sp)
# library(lattice) # required for trellis.par.set():
# #trellis.par.set(sp.theme()) # sets color ramp to bpy.colors()
# library(corrgram)
# #library(gdalUtils)
# library(corrgram)
library(dplyr)
# library(ape)
library(lavaan)
library(corrgram)
name <- function(x) { as.data.frame(names(x))} 

d <- read.csv("calib.data-4.2.csv")[,-1]
p <- unique(d)
p$hor <- as.factor(p$hor)
p <- p[p$top == 0 | p$hor == "B" | p$hor == "C", ]

# u <- d[,35:80]
# corrgram(u)
# names(u) <- gsub(x = names(u),pattern = "_mean",replacement = "m")
# names(u) <- gsub(x = names(u),pattern = "_std",replacement = "s")
# name(u)
name(p)
p$lstm <- p$lstm/100
p$lstsd <- p$lstsd/10
p$dem <- p$dem/100
p$river <- p$river/100000
p$XX1 <- p$XX1/1000
p$XX2 <- p$XX2/1000
p$XX3 <- p$XX3/1000
p$XD1 <- p$XD1/1000
p$XD2 <- p$XD2/1000
p$XD3 <- p$XD3/1000
p$X <- (p$X-5500000)/1000000
p$Y <- (p$Y-6000000)/1000000
p$wdist <- log(p$wdist)/10
p$wdist[p$wdist == -Inf] <- 0
p$clay <- p$clay/100
p$silt20 <- p$silt20/100
p$CEC <- p$CEC/10

name(p)
p1 <- p[,c(1,3,8,12:14,35:42)]
p1 <- p1[complete.cases(p1),]
var(p$dem[p$])
model1 <- '
# measurement model
cl =~ 1*clay
cec =~ 1*CEC
terrain =~ river + dem
evi =~ XX1 + XX2 + XX3
lu =~ XD1 + XD2 + XD3

# structural moel
terrain ~ X + Y
cl ~  X  + terrain + lu
cec ~  X + Y + cl + terrain + evi

# measurement error
river ~~ a*river
a > 0

# covariance 
XX2 ~~ XD2 + XX3
XX3 ~~ XD3 + XD2 
XX1 ~~   X
X  ~~ Y
dem ~~ Y 
clay ~~ X

# intercept
X ~ -0.025 * 1
Y ~ 0.225 * 1

# variance
X ~~ 0.0027 * X
Y ~~ 0.0014 * Y
'

model2 <- '
# measurement model
cl =~ 1*clay
cec =~ 1*CEC

# structural moel
#river + dem ~ X + Y
cl ~  X  + river + dem + XX1 + XX2 + XX3 
cec ~  X + Y + cl + dem + XX1 + XX2 + XX3 + XD1 + XD2 + XD3

# measurement error
# clay ~~ 0.1*clay
# CEC ~~ 0.1*CEC

# model variance
# cl ~~ a*cl
# a>0

# covariance 
XD1 ~~ XD3 + XX1
XD3 ~~ XX1
X ~~ XX1 + XD3
XD2 + XX3 ~~ XX2
XD1 ~~ XX3
XX1 ~~ XX3
XX3 ~~ XD3
dem ~~ XX1 + XD1 + XX3 
X + Y ~~ dem + river
river ~~ dem
dem ~~   XD3

# intercept
X ~ -0.025 * 1
Y ~ 0.225 * 1
XX1 ~ 5.414 * 1
XX2 ~ 2.706 * 1
XX3 ~ 3.303 * 1
XD1 ~ 1.258 * 1
XD2 ~ 0.780 * 1
XD3 ~ 1.210 * 1

# variance
X ~~ 0.0027 * X
Y ~~ 0.0014 * Y
XX1 ~~ 0.617 * XX1
XX2 ~~ 0.297 * XX2
XX3 ~~ 0.363 * XX3
XD1 ~~ 0.126 * XD1
XD2 ~~ 0.063 * XD2
XD3 ~~ 0.167 * XD3
'
# cor(p1[3:14])

fit1 <- sem(model2, data = p[p$hor=="B",], meanstructure = T,# group = "hor",group.label = c("A","B", "C"),
            fixed.x = F, verbose=TRUE)#, test = "bollen.stine") for low number of samples
inspect(fit1,"cov.lv")
inspect(fit1,"theta")
summary(fit1, fit.measures=F, standardized = F, rsquare = T)
varTable(fit1)
partable(fit1)

modindices(fit1)
modindices(fit1)[modindices(fit1)$op == "~~" & modindices(fit1)$mi >10,]

mean(p$Y[p$hor=="A"])

evi.mod <- '
evi =~ XX1 + XX2 + XX3 + XD2 + XD3 

'
fit.evi <- sem(evi.mod, data = p[p$hor=="B",], meanstructure = T,# group = "hor",group.label = c("A","B", "C"),
            fixed.x = F, verbose=TRUE)
