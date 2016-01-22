setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration")
rm(list=ls())
#install.packages('reshape')
# install.packages("lavaan", repos="http://www.da.ugent.be", type="source")
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
library(corrgram)
library(dplyr)
library(ape)
name <- function(x) { as.data.frame(names(x))} 

d <- read.csv("calib.data-4.2.csv")[,-1]

#Let's look at the spatial structure
library(ggplot2)

qplot(X, Y, data=d[d$hor == "B",], size=clay, color=river, alpha= 0.5) +
  theme_bw(base_size=18) + 
  scale_size_continuous("Clay %", range=c(0,10)) + 
  scale_color_gradient("Parana river", low="#0A0A2A", high="#F2F5A9")

######## Argiluviation ##########

a <- group_by(d, id.p)
a0 <- a[a$top==0,]
ab1 <- summarise(a[a$hor=="B",],max(clay, na.rm = T))
ab0 <- merge(ab1, a0, by = "id.p", all.x = T)
names(ab0)[2] <- "clay.B"

ab0$clay.ratio <- ab0$clay.B / ab0$clay

corrgram(ab0[,c(21:35,87)],lower.panel=panel.pie, upper.panel=panel.conf, diag.panel=panel.minmax)


step(lm(log(clay.ratio) ~ dem + river + wdist + maxc + mrvbf + slope + 
          twi + vdchn + water + lstm + lstsd + evim + evisd + 
#           X001_mean + X017_mean + X033_mean + X049_mean + X065_mean + 
#           X081_mean + X097_mean + X113_mean + X129_mean + X145_mean + 
#           X161_mean + X177_mean + X193_mean + X209_mean + X225_mean + 
#           X241_mean + X257_mean + X273_mean + X289_mean + X305_mean + 
#           X321_mean + X337_mean + X353_mean +
#           X001_std + X017_std + X033_std + X049_std + X065_std + X081_std + X097_std +
#           X113_std + X129_std + X145_std + X161_std + X177_std +
#           X193_std + X209_std + X225_std + X241_std + X257_std +
#           X273_std + X289_std + X305_std + X321_std + X337_std +
          X353_std + X + Y +
          XX1 + XX2 + XX3, ab0),direction = "both")

summary(lm(formula = log(clay.ratio) ~ dem + river + lstm + lstsd + evisd, 
           data = ab0))

summary(lm(formula = clay.B ~ dem + river + wdist + slope + lstm + lstsd + 
             XX2, data = ab0))

summary(lm(formula = clay ~ river + wdist + evisd, data = ab0))


ab0$log.clay.ratio <- log(ab0$clay.ratio)*100
ab0$river <- ab0$river/10000
ab0$dem <- ab0$dem/10
ab0$lstm <- ab0$lstm-273

ab1 <- ab0[,c(1,2,22,23,24,30,32:35,82:87)]
ab1 <- ab1[complete.cases(ab1),]
corrgram(ab1,lower.panel=panel.pie, upper.panel=panel.conf, diag.panel=panel.minmax)
hist(log(ab0$clay.ratio))

# summary(lm(clay.B~dem + river + lstm + lstsd, ab1))
# summary(lm(clay.ratio~dem + river + lstm + lstsd, ab1))

library(lavaan)
var(ab1$dem, na.rm = T)
var(ab1$river, na.rm = T)
var(ab1$lstm, na.rm = T)
var(ab1$lstsd, na.rm = T)
mean(ab1$dem, na.rm = T)
mean(ab1$river, na.rm = T)
mean(ab1$lstm, na.rm = T)
mean(ab1$lstsd, na.rm = T)
ab1$water[ab1$water>0] <- 1
ab1$water <- as.factor(ab1$water)

argil.mod <- '
argi =~  clay.B + clay.ratio 
argi ~ dem + river + lstm + lstsd 

# variances external drivers are fixed
dem ~~ 3.886*dem
river ~~ 14.373*river
lstm ~~ 0.183*lstm
lstsd ~~ 0.066*lstsd

# intercepts external drivers are fixed
dem ~ 6.020 * 1
river ~ 7.2659 * 1
lstm ~  23.778 * 1
lstsd ~ 7.719 * 1

dem ~~ river + lstsd + lstm
river ~~  lstm + lstsd
lstm ~~ lstsd
clay.B ~~ river
clay.ratio ~~ clay.ratio
clay.B ~~     5*clay.B 
#indirect := a*b
# clay.B ~~     clay.ratio 
'
fit1 <- sem(argil.mod, ab1, fixed.x = FALSE, meanstructure = T)
summary(fit1, rsquare = T, modindices = F, standardized = T)
modindices(fit1)
parameterEstimates(fit1)
#as.data.frame(fitMeasures(fit1, "all"))
# coef(fit1, "all") 

thing.to.save <- cbind(clay.B.est = rep(NA, nrow(ab1)),
                       argi = rep(NA, nrow(ab1)),
                       Xsqr = rep(NA, nrow(ab1)),
                       clay.B.mean.res = rep(NA, nrow(ab1)))  
for (i in 1:nrow(ab1)) {
  out <- sem(argil.mod, ab1[-i,], fixed.x = FALSE, meanstructure = T)
  #thing.to.save[i] <- inspect(object = out, "cov.lv")
  thing.to.save[i,1] <- parameterEstimates(out)$est[parameterEstimates(out)$label=="a"]
#   abx <- ab1[i,]
#   abx[,c(2,15)] <- c(10,1.5)
  thing.to.save[i,2] <- mean(lavPredict(object = out,type = "lv"))
  thing.to.save[i,3] <- fitMeasures(out, "chisq")
  thing.to.save[i,4] <- resid(out)$mean[2]
}
library(psych)
pairs.panels(thing.to.save,lm = T)
ab3 <- cbind(ab1,thing.to.save)
name(ab3)
pairs(ab3[,c(2,15:17)])^2
cor(ab3[,c(2,15:17)])^2
hist(thing.to.save[,2],30)
#ab1 <- ab1[,-6]
#cor(cbind(ab1,residuals(fit1, "casewise")))^2
# predict argi and residuals
ab2 <- cbind(ab1, argi = lavPredict(fit1,"lv"),
             clay.B.res = resid(fit1,"casewise")[,2],
             clay.ratio.res = resid(fit1,"casewise")[,1]
             )

pairs(ab2[,c(18,19)])
cor(ab2[,c(18,19)])^2

ab2$clay.B.pred <- ab2$clay.B + ab2$clay.B.res
name(ab2)
pairs(ab2[,c(20,2)])
cor(x = ab2$clay.B.pred, y = ab2$clay.B)^2

ab2$clay.ratio.pred <- ab2$clay.ratio + ab2$clay.ratio.res
name(ab2)
pairs(ab2[,c(21,16)])
cor(x = ab2$clay.ratio.pred, y = ab2$clay.ratio)^2

#Let's look at the spatial structure
library(ggplot2)

qplot(X, Y, data=ab2, color=clay.ratio.res,size=I(5)) +
  theme_bw(base_size=18) + 
  #scale_size_continuous("DEM", range=c(0,10)) + 
  scale_color_gradient("Argilluviation", low="#0A0A2A", high="#F2F5A9")

###### END Argiluviation #########


name(a)
b <- a[,c(1:34,83:85, 81, 82)]
corrgram(b[b$hor!="A",],lower.panel=panel.pie, upper.panel=panel.conf, diag.panel=panel.minmax)
b$resist <- (-1)*((b$resist/1000)^0.5)

mean(b$tb[b$hor!="A"], na.rm = T)
mean(b$CEC[b$hor!="A"], na.rm = T)
mean(b$phw[b$hor!="A"], na.rm = T)
mean(b$phkcl[b$hor!="A"], na.rm = T)
mean(b$resist[b$hor!="A"], na.rm = T)

summary(b[b$hor!="A",c(21,22,28,31:37)])
b$tb <- b$tb/10
b$dem <- b$dem/10
b$river <- b$river/10000
b$lstm <- b$lstm/100
b$evim <- b$evim/1000
b$evisd <- b$evisd/1000
b$XX1 <- b$XX1/10
b$XX2 <- b$XX2/1000
b$XX3 <- b$XX3/1000

mod2 <- '
pH =~ phw + phkcl 
#Cat =~ CEC + tb
Env =~ evisd + XX2
Ter =~ dem + river

pH ~ Env + Ter
#Cat ~ pH + Env + Ter
'

fit2 <- sem(mod2, b[b$hor!="A",], fixed.x = FALSE, meanstructure = T)
varTable(fit2)
summary(fit2, rsquare = T, modindices = F, standardized = T)
modindices(fit2)















