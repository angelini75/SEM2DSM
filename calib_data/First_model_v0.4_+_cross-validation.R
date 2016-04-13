rm(list=ls())
# install.packages("lavaan")
# install.packages("lavaan", repos="http://www.da.ugent.be", type="source")
name <- function(x) { as.data.frame(names(x))} # as.data.frame(.Primitive("names"))
library(lavaan)
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")
#original <- read.csv("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/OLD/calib.data-1.0.csv")
d <- read.csv("calib.data-2.1.csv")[,-1]

############### PRE-PROCESSING ################## 
name(d)
names(d)[c(5,6,9,10)] <- c("tb.A","sat.A", "oc.A","bt")
d <- d[!is.na(d$oc.A),]
d$sat.A[d$id.p==480] <- 88 #error in dataset
# statistics of calibration data
D<-d[,4:10]
D$esp.A <- log10(D$esp.A)
D$esp.B <- log10(D$esp.B)
D <- D[,c(1,6,2:5,7)]
d.stat<- matrix(data = NA,nrow = 6,ncol = 7,
                dimnames = list(c("Min","Median","Mean", "Max", "SD","SS"),names(D)))
d.stat[1,]<- apply(X = D,FUN = min,2, na.rm=T) # 2 means by column
d.stat[2,]<- apply(X = D,FUN = median,2, na.rm=T)
d.stat[3,]<- apply(X = D,FUN = mean,2, na.rm=T)
d.stat[4,]<- apply(X = D,FUN = max,2, na.rm=T)
d.stat[5,]<- apply(X = D,FUN = sd,2, na.rm=T)
for(i in 1:7){
  d.stat[6,i] <- sum((mean(D[,i], na.rm=T)-D[,i]) ^ 2, na.rm=T)
}
write.csv(d.stat, "d.stat.cal.csv")
stat.desc(D)


### assumptions
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
name(d)
#order of soil properties: thick, oc, tb, sat, esp.a, esp.b, bt
d <- d[,c(1:3,4,9,5:8,10,11:29)]

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

#### Third Model####
third_model <- '
# measurement model
thick.Ar =~ 1*thick.A
oc.Ar =~ 1*oc.A
tb.Ar =~ 1*tb.A
sat.Ar =~ 1*sat.A
esp.Ar =~ 1*esp.A
esp.Br =~ 1*esp.B
btr =~ 1*bt

# structural model
thick.Ar ~  dem + wdist + mrvbf + vdchn + twi + river + slope + maxc + evim + evisd

oc.Ar ~     lstm +  lstsd + evim + evisd + dem + wdist + mrvbf + vdchn + twi +
            esp.Br + esp.Ar + btr + thick.Ar
tb.Ar ~     evim + evisd + lstm + lstsd + dem + wdist + mrvbf + vdchn + twi + river + 
            oc.Ar + btr
sat.Ar ~    evim + evisd + lstm + lstsd + dem + wdist + mrvbf + vdchn +  twi + river + 
            tb.Ar + oc.Ar                            
esp.Br ~    lstm +  lstsd + dem + wdist + mrvbf + vdchn + twi + river

esp.Ar ~    lstm +  lstsd + dem + wdist + mrvbf + vdchn + twi + river + 
            esp.Br 
btr ~       lstm +  lstsd + wdist + vdchn + twi + dem + river + mrvbf +
            esp.Br + esp.Ar

# measurement error
thick.A ~~  0.25*thick.A
oc.A ~~     0.20*oc.A
tb.A ~~     0.20*tb.A
sat.A ~~    0.20*sat.A
esp.A ~~    0.20*esp.A
esp.B ~~    0.10*esp.B
bt  ~~      0.25*bt
'
##### fitting ####
#fit3<- sem(third_model, d,meanstructure = T,std.lv = T, ordered = c("is.E","is.caco3","is.hydro"))

library(utils)
pred <- as.data.frame(NA)
pred[,1:38] <- NA
names(pred)[1:29] <- names(d)
names(pred)[32:38] <- c("thick.Ar","oc.Ar","tb.Ar","sat.Ar","esp.Ar","esp.Br","btr") #names of soil properties

pb = txtProgressBar(min = 0, max = length(d[,1]), initial = 0, style = 3)

Var <- matrix(nrow = 0, ncol = 7, dimnames = list(NULL,
              c("thick.Ar","oc.Ar","tb.Ar","sat.Ar","esp.Ar","esp.Br","btr")))

resids <- pred[32:38]
theta <- pred[,c(32:38)]
a <- pred[,c(32:38)] #observed
b <- pred[,c(32:38)] #predicted
v <- pred[,c(32:38)] #variance(s)
resids <- pred[,c(32:38)] #residuals

for (i in 1:length(d[,1])) {
calib <- d[-i,]
pred[i,] <- d[i,]
fit3 <- lavaan(model = third_model, data = calib,
               model.type = "sem", meanstructure = "default",
               int.ov.free = FALSE, int.lv.free = FALSE, fixed.x = "default",
               orthogonal = FALSE, std.lv = FALSE, 
               parameterization = "default", auto.fix.first = F,
               auto.fix.single = T, auto.var = T, auto.cov.lv.x = FALSE,
               auto.cov.y = T, auto.th = T, auto.delta = FALSE,
               std.ov = FALSE, missing = "default",
               constraints = "", estimator = "ML",
               zero.cell.warn = TRUE, start = "default")
#summary(fit3)
# modi<-summary(fit3, standardized=F, modindices = T, fit.measures=F) 
# modi<-modi[modi$mi>3 & !is.na(modi$mi),]
# modi

################# Matrices ##############
##setting up matrices
B <- inspect(fit3, "est")$beta[1:7,1:7] # matrix of coeff. latent state variables
I <- diag(nrow = 7, ncol = 7) # Identity matrix
A <- inspect(fit3, "est")$beta[1:7,8:19] # matrix of coeff of external drivers
V <- inspect(fit3, "est")$psi[1:7,1:7] # matrix of predicted error variance
Th <- inspect(fit3, "est")$theta[1:7,1:7] # matrix of measurement error
IB <- solve(I - B)
################# Running Prediction ###########################
# (IB%*%A%*%p) product of matrices per pixel (equation 4 paper)
  p = as.vector(as.matrix(pred[i,c(15,17,19,22,21,16,20,18,26,27,24,25)]))
  p = matrix(p, nrow = 12, ncol = 1)
  pred[i,32:38] = t(IB %*% A %*% p) # key equation
  # calculate standarised squared standard error
  ## theta is standarised squared standard error
  ## theta = ((observed-predicted)²)/error variance=(standard error^2)
  a[i,] <- pred[i,c(4:10)] # observed values
  b[i,] <- pred[i,c(32:38)] # predicted values
  v[i,] <- diag(IB%*%V%*%t(IB)+Th) # error variance
  resids[i,] <- a[i,] - b[i,] # residuals
  theta[i,] <- (resids[i,]^2)/v[i,] # theta 

##### Estimation confidence interval ## No stapial #######
Var.n <- IB %*% V %*% t(IB) # diagonal vaues are variance error
Var <- rbind(Var, diag(Var.n))
#bar time
setTxtProgressBar(pb, i)
}
summary(Var)
Var <- apply(Var, MARGIN = 2, FUN = mean)

name(pred)
res <- pred[,c(4,4+28,5,5+28,6,6+28,7,7+28,8,8+28,9,9+28,10,10+28)]
names(res) <- rep(names(pred)[4:10], each=2)
names(res)[c(2,4,6,8,10,12,14)] <- paste(names(res)[c(2,4,6,8,10,12,14)], "p", sep = ".")

M <- N[rep(1:7, each=2),]

# Un-standardize
for (i in 1:14) { 
  res[,i] <- res[,i] * M[i,3] + M[i,2]
}
#back-transform
# res$esp.A <- 10 ^ (res$esp.A)
# res$esp.A.p <- 10 ^(res$esp.A.p + (Var[5] ^ 2) * 0.5)
# res$esp.B <- 10 ^ (res$esp.B)
# res$esp.B.p <- 10 ^ (res$esp.B.p + (Var[6] ^ 2) * 0.5)

#par()
par(mfrow = c(3, 3), pty="s",mai=rep(0.7,4))
l <- c(2,4,6,8,10,12,14)
for (i in l) {
lim = c(min(res[,i - 1]), max(res[,i - 1]))
plot(res[,i] ~ res[,i - 1], main = paste(names(res)[i - 1]), xlab = "measured",
     ylab = "predicted", col = "dark red", xlim = lim, ylim = lim)
abline(0,1)
abline(lm(res[,i] ~ res[,i - 1]), col = "blue")
}


report <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, SS = NA, mean_theta = NA, median_th = NA)
for (i in l) {
######################## ME <- mean error 
  ME  <-  mean(res[,i] - res[,i - 1])
############################################ RMSE (root mean squared error)
  RMSE <- sqrt(mean((res[,i] - res[,i - 1]) ^ 2))
  MSE <- mean((res[,i] - res[,i - 1]) ^ 2)
############################################ SS (Sum of squares)
  SS <- sum((res[,i-1] - res[,i]) ^ 2)
# fill report table
  report[i / 2,1:4] <- c(names(res)[i - 1], ME, RMSE, SS)
}


############################################ theta mean and median 
for(i in 1:7){
  report$mean_theta[i] <- mean(theta[,i])
  report$median_th[i] <- median(theta[,i])
}
report
#d.stat <- read.csv("summary.calibdata.csv")
report$R2 <- 1 - (as.numeric(report$SS) / N$SStot[1:7])
report


write.csv(report, "/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/report.SEM.cros-val.csv")

# #########################################################################################################################
# #########Plot
# library(semPlot)
# #matrix to arrange nodes at graph
# layout <-as.matrix(read.csv("matrix_semplot.csv")[,-1])
# # plot sem using layout=layout
# 
# layout[]
# 
# semPaths(fit3,what = "std",whatLabels = "no", layout = layout,sizeLat = 4, cut =0.35,
#          sizeInt2 = 2,sizeMan =3.5, sizeLat2 = 2, sizeInt = 1,sizeMan2 = 1.5,nCharNodes=0, font=3,
#          edge.width = 2,esize=1.5, asize=1,intercepts = F, reorder = F,equalizeManifests =T, residuals = F,layoutSplit=F,
#          structural = F, exoCov = F, exoVar=F,cardinal = F,style = "lisrel",#color = c("orange","blue"),
#          manifests = c("tb.A", "sat.A", "evim", "evisd", "lstm", "lstsd", "dem","bt", "wdist", "mrvbf", "vdchn", 
#                        "twi", "river", "thick.A", "slope", "maxc", "oc.A", "esp.B", "esp.A"))
# 
# semPlotModel(fit3)
# semPaths(fit3)
# 
# str(fit3)
# 
# write.csv(report, "/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/report_cross-validation.csv")
# 
# 
# semPaths




