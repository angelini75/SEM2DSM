
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration")
rm(list=ls())
#install.packages('reshape')
#install.packages("lavaan")#, repos="http://www.da.ugent.be", type="source")
# install.packages("semTools")
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
# install.packages('lavaan.survey')
# library(lavaan.survey)
library(dplyr)
# library(ape)
library(lavaan)
library(semTools)
library(corrgram)
name <- function(x) { as.data.frame(names(x))} 
####################


d <- read.csv("calib.data-4.2.csv")[,-1]
p <- unique(d)
p$hor <- as.factor(p$hor)
p <- p[p$top == 0 | p$hor == "A" | p$hor == "B" | p$hor == "C", ]
p$hor[p$top==0] <- "A"
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
p$vdchn <- p2$vdchn/10

name(p)
p1 <- p[,c(1,3,8,12:14,35:42)]
p1 <- p1[complete.cases(p1),]

skew(p$clay[complete.cases(p$clay)])
skew(p$OC[complete.cases(p$OC)])
skew(p$CEC[complete.cases(p$CEC)])
p2 <- indProd(data = p,var1 = 21:22,var2 = 40:42, match = F, meanC=FALSE, 
                    residualC=TRUE, doubleMC=FALSE)
p3 <- merge(x=p2[p2$hor=="A",1:15],y = p2[p2$hor=="B",1:15], by= "id.p", all = TRUE )
p3 <- merge(x=p3,y = p2[p2$hor=="C",1:15], by= "id.p", all = TRUE )
names(p3) <- gsub(pattern = ".x$",replacement = ".A",x = names(p3))
names(p3) <- gsub(pattern = ".y$",replacement = ".B",x = names(p3))
names(p3)[30:43] <- paste0(names(p3)[30:43],".C")

p3 <- merge(x = p3, y = unique(p2[,c(1,21:42)]), by = "id.p", all = TRUE)
names(p3)[41]<- "clay.C"

######
model1 <- '
# measurement model
cl =~ 1*clay
cec =~ 1*CEC
terrain =~ river + dem
evi =~ XX1 + XX2 + XX3
lu =~ XD1 + XD2 + XD3

# structural moel
terrain ~ X + Y
cl ~  X  + terrain + 0.6184101*lu
cec ~  X + Y + cl + terrain + 0.1466529*evi

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
######
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
######
model3 <- '
river ~ X + Y 
dem ~ river + X
river~~dem 
#X ~~   Y

lstm ~ river + X
lstsd ~ lstm + X + Y + dem

clay ~ dem + river + lstsd + X + lstm
CEC ~  dem + river + lstsd + X + clay
'

fit1 <- sem(model3, data = p2[p2$hor=="B",], meanstructure = T,# group = "hor",group.label = c("A","B", "C"),
            fixed.x = T, verbose=F)#, test = "bollen.stine") for low number of samples
inspect(fit2,"cov.lv")
inspect(fit2,"theta")
summary(fit1, fit.measures=F, standardized = F, rsquare = T)

# X ~~  -0.0005693461*   Y
# X ~~ 0.002739*X
# Y ~~ 0.00138*Y
# X ~ -0.020*1
# Y ~ 0.226*1
#######
model4 <- '
cl =~ 1*clay
cec =~ 1*CEC
oc =~ 1*OC
oc ~ cl
cl ~ vdchn + water + lstsd
cec ~  cl + vdchn + lstsd
cec ~~ 0*oc
'
fit2 <- sem(model4, data = p2[p2$hor=="C",], meanstructure = T,# group = "hor",group.label = c("A","B", "C"),
            fixed.x = T, verbose=F)#, test = "bollen.stine") for low number of samples
summary(fit2, fit.measures=F, standardized = F, rsquare = T)
a<-modindices(fit2)[modindices(fit2)$mi>0.6,]


varTable(fit2)
partable(fit2)
#####
mod.groups <- '
A:
cl =~ 1*clay
cec =~ 1*CEC
river ~ a1*X + a2*Y + a3*dem
dem ~ a4*Y + a5*X
lstm ~ a8*river + a9*X 
lstsd ~ a11*lstm + a12*X + a13*Y + a14*dem + vdchn


cl ~ dem + lstm + lstsd + X + Y + river + vdchn
cec ~  dem + X + Y + cl

B:
cl =~ 1*clay
cec =~ 1*CEC
river ~ a1*X + a2*Y + a3*dem
dem ~ a4*Y + a5*X
lstm ~ a8*river + a9*X 
lstsd ~ a11*lstm + a12*X + a13*Y + a14*dem 
river ~~ a16*vdchn
cl ~1

cl ~ dem + river + lstsd + X + lstm
cec ~  dem + river + lstsd + X + clay 

C:
cl =~ 1*clay
cec =~ 1*CEC
river ~ a1*X + a2*Y + a3*dem + vdchn
dem ~ a4*Y + a5*X
lstm ~ a8*river + a9*X + Y
lstsd ~ a11*lstm + a12*X + a13*Y + a14*dem + vdchn

cl ~1

cl ~ dem + lstm + lstsd + X + Y + river + vdchn
cec ~  dem + X + Y + cl 
'

fit2 <- sem(mod.groups, data = p2, meanstructure = T, group = "hor",group.label = c("A","B", "C"),
            fixed.x = T, verbose=F)#, test = "bollen.stine") for low number of samples
summary(fit2, fit.measures=T, standardized = F, rsquare = T)
modindices(fit2)[modindices(fit2)$mi>3000000,]
findRM
#####
mod.gr <- '
A:
cl =~ 1*clay
cec =~ 1*CEC
oc =~ 1*OC
oc ~ evisd + evim + cl + X + Y + bottom + dem + mrvbf 
cl ~ dem + X + river + evisd + bottom
cec ~  Y + cl + wdist + water + lstm + lstsd + oc
clay ~~ clay
#cl ~~    oc
#CEC ~~ CEC

B:
cl =~ 1*clay
cec =~ 1*CEC
oc =~ 1*OC
oc ~ evisd + lstm + cl + X + Y + bottom + dem + vdchn
cl ~ dem + vdchn + X + river + evisd + bottom + lstm + lstsd
cec ~  X + Y + cl + evim + lstm + lstsd + oc + bottom

C:
cl =~ 1*clay
cec =~ 1*CEC
oc =~ 1*OC
oc ~ lstm + cl
cl ~ dem + X + Y + bottom + lstm + lstsd + vdchn
cec ~  X + Y + cl + dem + lstm + bottom + vdchn
'
fit3 <- sem(mod.gr, data = p2, meanstructure = T, group = "hor",group.label = c("A","B", "C"),
            fixed.x = T, verbose=F)#, test = "bollen.stine")# for low number of samples
summary(fit3, fit.measures=T, standardized = F, rsquare = T)
modindices(fit3)[modindices(fit3)$mi>0,]

############

mod.p3 <- '
cl =~ clay.A + clay.B + clay.C
oc =~ OC.A  
cec =~ CEC.A + CEC.B + CEC.C

oc ~~ a*oc
OC.A ~~ a*OC.A
a>0

cl ~ X + Y + river + dem
oc ~ X + Y + river + dem + cl
cec ~ oc + cl 

clay.A ~~ CEC.A
OC.A ~~ CEC.A
clay.B ~~ CEC.B
clay.C ~~ CEC.C

clay.A + clay.B ~~ clay.C
clay.A ~~ clay.B
CEC.A ~~ CEC.C
'
fit3 <- sem(mod.p3, data = p3, meanstructure = T,# group = "hor",group.label = c("A","B", "C"),
            fixed.x = T, verbose=F)#, test = "bollen.stine")# for low number of samples
summary(fit3, fit.measures=T, standardized = F, rsquare = T)
modindices(fit3)[modindices(fit3)$mi>3,]

inspect(fit3,"cov.lv")
inspect(fit3,"theta")


























