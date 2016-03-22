rm(list=ls())
# install.packages("lavaan")
# install.packages("lavaan", repos="http://www.da.ugent.be", type="source")
name <- function(x) { as.data.frame(names(x))} # as.data.frame(.Primitive("names"))
library(lavaan)
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")
#original <- read.csv("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/OLD/calib.data-1.0.csv")
d <- read.csv("calib.data-2.1.csv")[,-1]



############### PRE-PROCESSING ################## 
names(d)
names(d)[c(5,6,9,10)] <- c("tb.A","sat.A", "oc.A","bt")
d$sat.A[d$id.p==480] <- 88 #error in dataset
# statistics of calibration data
D<-d[,4:10]
d.stat<- matrix(data = NA,nrow = 7,ncol = 7,
                dimnames = list(c("Min","Median","Mean", "Max", "SD","SS","N"),names(D)))
summary(D)
library(Hmisc)
describe(D, exclude.missing=TRUE,digits=4)
#install.packages("pastecs")
library(pastecs)
d.stat <- stat.desc(D,basic = T)

write.csv(d.stat, "/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/d.stat.cal.csv")
####



d <- d[!is.na(d$oc.A),]
d$sat.A[is.na(d$sat.A)] <- 100 # it is assumed 100% saturation when CaCO3 is present 
d$d.caco3[is.na(d$d.caco3)] <- 300 # it is assumed that CaCO3 is very deep when it is absent within the solum 
summary(d$tb.A[d$is.caco3 == 1], omit.na = TRUE)
d$tb.A[is.na(d$tb.A)] <- 20

# plot(d$esp.A[d$is.caco3==1]~d$d.caco3[d$is.caco3==1], omit.na=T)
# boxplot(d$d.caco3[is.na(d$esp.A)])
summary(d$esp.A[d$d.caco3 < 15], omit.na = TRUE)
d$esp.A[is.na(d$esp.A)] <- 16.3

# plot(d$esp.B[d$is.caco3==1]~d$esp.A[d$is.caco3==1])
summary(d$d.caco3[is.na(d$esp.B)])
summary(d$esp.B[d$d.caco3 < 30], omit.na = TRUE)
d$esp.B[is.na(d$esp.B)] <- 30.9

nas<-d[!complete.cases(d),]


# transformation
d$esp.A <- log10(d$esp.A)
d$esp.B <- log10(d$esp.B)


# save mean and sd
N <- data.frame(mean = rep(0,20), sd = rep(0,20), SStot=rep(0,20))
dm <- d[,c(4:10,15:27)]
for(i in 1:20){
  N$mean[i] <- mean(dm[,i])
  N$sd[i] <- sd(dm[,i])
  N$SStot[i] <- sum((mean(dm[,i])-dm[,i])^2)
  rownames(N)[i] <- names(dm)[i]
}
N <- data.frame(name=rownames(N),mean=N$mean,sd=N$sd, SStot=N$SStot)

# normalization
n <- c(4:10,15:27)
for(i in n){
  d[,i] <- (d[,i] - mean(d[,i])) / sd(d[,i])
}

############### FITTING MODEL ######################
name(d)
d <- d[,c(1:4,9,5,6:8,10,11:29)]
library(utils)
pred <- as.data.frame(NA)
pred[,1:36] <- NA
names(pred)[1:29] <- names(d)
names(pred)[30:36] <- c("thick.Ar","oc.Ar","tb.Ar","sat.Ar","esp.Ar","esp.Br","btr") #names of soil properties

theta <- pred[,c(1,30:36)]
resids <- pred[,c(1,30:36)]

for (i in 1:length(d[,1])) {
calib <- d[-i,]
pred[i,] <- d[i,]
tk <- lm(formula = thick.A ~ dem + evisd + maxc, data = calib)
oc <- lm(formula = oc.A ~ evisd + lstm + lstsd + mrvbf, data = calib)
tb <- lm(formula = tb.A ~ evim + evisd + lstm + lstsd + maxc + river + wdist, data = calib)
sa <- lm(formula = sat.A ~ dem + evim + evisd + maxc + river, data = calib)
ea <- lm(formula = esp.A ~ evisd + lstsd + mrvbf + twi + vdchn, data = calib)
eb <- lm(formula = esp.B ~ evisd + lstsd + mrvbf + river + twi + vdchn, data = calib)
bt <- lm(formula = bt ~ evisd + maxc + river + vdchn, data = calib)
  

models <- list(tk,oc,tb,sa,ea,eb,bt)




################# Running Prediction ###########################
## estimate residuals
## theta is standarised squared standard error
## theta is = ((observed-predicted)²)/variance error=(standard error^2)
 
for(j in 1:7){
  pred[i,29+j] <- predict(models[[j]],pred[i,])
  resids[i,1+j] <- pred[i,29+j] - d[i,3+j]
  theta[i,1+j] <- ((resids[i,1+j])^2)/summary(tk)$sigma^2
  }
 }

name(pred)
res <- pred[,c(1,4:10,30:36)]
M <- N[c(1,6,2:5,7),]
# Un-standardize
for (i in 2:8) { 
  res[,i] <- res[,i] * M$sd[i-1] + M$mean[i-1]
}
for (i in 9:15) { 
  res[,i] <- res[,i] * M$sd[i-8] + M$mean[i-8]
}

#back-transform
# res$esp.A  <- 10 ^ (res$esp.A)
# res$esp.Ar <- 10 ^ (res$esp.Ar + var(res$esp.Ar) * 0.5)
# res$esp.B  <- 10 ^ (res$esp.B)
# res$esp.Br <- 10 ^ (res$esp.Br + var(res$esp.Br) * 0.5)

#par()
par(mfrow = c(3, 3), pty="s",mai=rep(0.7,4))
l <- c(2:8)
for (i in l) {
lim = c(min(res[,i]), max(res[,i]))
plot(res[,i+7] ~ res[,i], main = paste(names(res)[i]), xlab = "measured",
     ylab = "predicted", col = "dark red", xlim = lim, ylim = lim)
abline(0,1)
abline(lm(res[,i+7] ~ res[,i]), col = "blue")
}


report<- data.frame(Soil_property = NA, ME = NA, RMSE = NA, SS = NA, mean_theta = NA, median_th = NA)


for (i in 2:8) {
######################## ME <- mean error 
  ME <- mean(res[,i] - res[,i + 7])
  RMSE <- sqrt(mean((res[,i] - res[,i + 7]) ^ 2))
  MSE <- mean((res[,i] - res[,i + 7]) ^ 2)
  SS <- sum((res[,i] - res[,i + 7]) ^ 2)
  mean_th <- mean(theta[,i])
  median_th <- median(theta[,i])
# fill report table
  report[i-1,1:6] <- c(names(res)[i], ME, RMSE, SS, mean_th, median_th)
}
report
#d.stat <- read.csv("summary.calibdata.csv")

# d.stat <- d.stat[,c(1,6,2:5,7)]
# D <- D[,c(1,6,2:5,7)]
M
for(i in 1:7){
report$R2[i] <- 1 - (as.numeric(report$SS[i]) / M$SStot[i])
}

report

write.csv(report, "/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/report.lm.cros-val.csv")

