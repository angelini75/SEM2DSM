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

# statistics of calibration data
D<-d[,4:10]
d.stat<- matrix(data = NA,nrow = 6,ncol = 7,
                dimnames = list(c("Min","Median","Mean", "Max", "SD","SS"),names(D)))
d.stat[1,]<- apply(X = D,FUN = min,2) # 2 means by column
d.stat[2,]<- apply(X = D,FUN = median,2)
d.stat[3,]<- apply(X = D,FUN = mean,2)
d.stat[4,]<- apply(X = D,FUN = max,2)
d.stat[5,]<- apply(X = D,FUN = sd,2)
for(i in 1:7){
  d.stat[6,i] <- sum((mean(D[,i])-D[,i]) ^ 2)
}

# transformation
d$esp.A <- log10(d$esp.A)
d$esp.B <- log10(d$esp.B)
d$is.hydro <- ordered(d$is.hydro)
d$is.E <- ordered(d$is.E)
d$is.caco3 <- ordered(d$is.caco3)

# hist(10^(d$esp.A),col = "lightblue")
# summary(d$thick.A)
##@## data normalization
N <- data.frame(mean = rep(0,20), sd = rep(0,20))

# save mean and sd

dm <- d[,c(4:10,15:27)]
for(i in 1:20){
  N$mean[i] <- mean(dm[,i])
  N$sd[i] <- sd(dm[,i])
  rownames(N)[i] <- names(dm)[i]
}

# normalization
n <- c(4:10,15:27)
for(i in n){
  d[,i] <- (d[,i] - mean(d[,i])) / sd(d[,i])
}

#step(lm(sat.A ~ dem+wdist+maxc+mrvbf+slope+twi+vdchn+lstm+lstsd+evim+evisd+river, d),direction = "both")
# summary(lm(formula = bt ~ dem + maxc + slope + lstm + evisd, data = d))
# summary(lm(formula = bt ~ dem + maxc + slope + lstm + evisd + tb.A + d.caco3 + river, data = d))
# oc.fit<- lm(formula = oc.A ~ lstm +  lstsd + evim + evisd + dem + wdist + mrvbf + vdchn + twi, data = d)
# summary(oc.fit)
# x<-predict(oc.fit,data=d)
# summary(lm(formula = tb.A ~ wdist + lstm + evim + evisd + river, data = d))
# summary(lm(formula = thick.A ~ dem + maxc + evisd, data = d))
# summary(lm(formula = esp.B ~ mrvbf + twi + vdchn + lstsd + evim + evisd + 
#              river, data = d))
# summary(lm(formula = esp.A ~ dem + mrvbf + twi + vdchn + lstsd + evim + 
#              esp.B, data = d))
# summary(lm(formula = sat.A ~ maxc + lstm + lstsd + evim + evisd, data = d))
# summary(lm(formula = d.caco3 ~ dem + mrvbf + vdchn + lstm + lstsd + evim + 
#             evisd + river, data = d))
N <- read.csv("N.csv")

############### FITTING MODEL ######################

#### Third Model####
third_model <- '
# measurement model
tb.Ar =~ 1*tb.A
sat.Ar =~ 1*sat.A
btr =~ 1*bt
oc.Ar =~ 1*oc.A
thick.Ar =~ 1*thick.A
esp.Br =~ 1*esp.B
esp.Ar =~ 1*esp.A

# structural model
tb.Ar ~     evim + evisd + lstm + lstsd + dem + wdist + mrvbf + vdchn + twi + river + 
            oc.Ar + btr
sat.Ar ~    evim + evisd + lstm + lstsd + dem + wdist + mrvbf + vdchn +  twi + river + 
            tb.Ar + oc.Ar                            
btr ~       lstm +  lstsd + wdist + vdchn + twi + dem + river + mrvbf +
            b*esp.Br + b*esp.Ar
oc.Ar ~     lstm +  lstsd + evim + evisd + dem + wdist + mrvbf + vdchn + twi +
            a*esp.Br + a*esp.Ar + btr + thick.Ar
esp.Ar ~    lstm +  lstsd + dem + wdist + mrvbf + vdchn + twi + river + evisd +
            esp.Br 
esp.Br ~    lstm +  lstsd + dem + wdist + mrvbf + vdchn + twi + river + evisd
            
thick.Ar ~  dem + wdist + mrvbf + vdchn + twi + river + slope + maxc + evim + evisd

#a == -0.5

# measurement error
thick.A ~~  0.25*thick.A
tb.A ~~     0.20*tb.A
sat.A ~~    0.20*sat.A
bt  ~~      0.25*bt
oc.A ~~     0.20*oc.A
esp.B ~~    0.10*esp.B
esp.A ~~    0.20*esp.A
# thick.A ~~  0.25*thick.A
# tb.A ~~     0.1*tb.A
# sat.A ~~    0.1*sat.A
# bt  ~~      0.1*bt
# oc.A ~~     0.1*oc.A
# esp.B ~~    0.10*esp.B
# esp.A ~~    0.1*esp.A
# intercepts
# tb.Ar ~1
'
##### fitting ####
#fit3<- sem(third_model, d,meanstructure = T,std.lv = T, ordered = c("is.E","is.caco3","is.hydro"))

library(utils)
pred <- as.data.frame(NA)
pred[,1:38] <- NA
names(pred)[1:29] <- names(d)
names(pred)[32:38] <- c("tb.Ar","sat.Ar","btr","oc.Ar","thick.Ar","esp.Br","esp.Ar") #names of soil properties

pb = txtProgressBar(min = 0, max = length(d[,1]), initial = 0, style = 3)

Var <- matrix(nrow = 0, ncol = 7, dimnames = list(NULL,
              c("tb.Ar","sat.Ar","btr","oc.Ar","thick.Ar","esp.Br","esp.Ar")))
for (i in 1:length(d[,1])) {
calib <- d[-i,]
pred[i,] <- d[i,]
fit3 <- lavaan(model = third_model, data = d,
               model.type = "sem", meanstructure = "default",
               int.ov.free = FALSE, int.lv.free = FALSE, fixed.x = "default",
               orthogonal = FALSE, std.lv = FALSE, 
               parameterization = "default", auto.fix.first = F,
               auto.fix.single = T, auto.var = T, auto.cov.lv.x = FALSE,
               auto.cov.y = T, auto.th = T, auto.delta = FALSE,
               std.ov = FALSE, missing = "default",
               constraints = "", estimator = "ML",
               zero.cell.warn = TRUE, start = "default")
summary(fit3)
modi<-summary(fit3, standardized=F, modindices = T, fit.measures=F) 
modi<-modi[modi$mi>3 & !is.na(modi$mi),]
modi

################# Matrices ##############
##setting up matrices
B <- inspect(fit3, "est")$beta[1:7,1:7] #matrix of coeff. latent state variables
I <- diag(nrow = 7, ncol = 7) #Identity matrix
A <- inspect(fit3, "est")$beta[1:7,8:19] #matrix of coeff of external drivers
V <- inspect(fit3, "est")$psi[1:7,1:7] #matrix of predicted error variance
IB <- solve(I - B)

################# Running Prediction ###########################

# (IB%*%A%*%p) product of matrices per pixel (equation 4 paper)
for (j in i) {
  p = as.vector(as.matrix(pred[j,c(26,27,24,25,15,17,19,22,21,16,20,18)]))
  p = matrix(p, nrow = 12, ncol = 1)
  pred[j,32:38] = t(IB %*% A %*% p) # key equation
  }
##### Estimation confidence interval ## No stapial #######
Var.n <- IB %*% V %*% t(IB) # diagonal vaues are variance error
Var <- rbind(Var, diag(Var.n))
#bar time
setTxtProgressBar(pb, i)
}
summary(Var)
Var <- apply(Var, MARGIN = 2, FUN = mean)

name(pred)
res <- as.data.frame(cbind(pred[,4], pred[,36], pred[,5], pred[,32], pred[,6], pred[,33],
         pred[,7], pred[,37], pred[,8], pred[,38], pred[,9], pred[,35], pred[,10], pred[,34]))
names(res) <- rep(names(pred)[4:10], each=2)
names(res)[c(2,4,6,8,10,12,14)] <- paste(names(res)[c(2,4,6,8,10,12,14)], "p", sep = ".")

M <- N[rep(1:7, each=2),]

# Un-standardize
for (i in 1:14) { 
  res[,i] <- res[,i] * M[i,3] + M[i,2]
}
#back-transform
res[,7] <- 10 ^ (res[,7])
res[,8] <- 10 ^(res[,8] + (Var[7] * N$sd[4] ^ 2) * 0.5)
res[,9] <- 10 ^ (res[,9])
res[,10] <- 10 ^ (res[,10] + (Var[6] * N$sd[5] ^ 2) * 0.5)

par()
par(mfrow = c(3, 3), pty="s",mai=rep(0.7,4))
l <- c(2,4,6,8,10,12,14)
for (i in l) {
lim = c(min(res[,i - 1]), max(res[,i - 1]))
plot(res[,i] ~ res[,i - 1], main = paste(names(res)[i - 1]), xlab = "measured",
     ylab = "predicted", col = "dark red", xlim = lim, ylim = lim)
abline(0,1)
abline(lm(res[,i] ~ res[,i - 1]), col = "blue")
}


report<- data.frame(Soil_property = NA, ME = NA, RMSE = NA, SS = NA)
N  <- as.numeric(length(res$thick.A))

for (i in l) {
######################## ME <- mean error 
  z = mean(res[,i] - res[,i - 1])
# variance of z
  Vz <- var(res[,i] - res[,i - 1])
# 95% confidence
# lowwer limit
  ll <- z - qt(0.975, N - 1) * sqrt(Vz)
# upper limit
  ul <- z + qt(0.975, N - 1) * sqrt(Vz)

  ME <- paste(round(ll, 3), "<", round(z, 3), "<", round(ul, 3))

############################################ RMSE (root mean squared error)
  RMSE <- sqrt(mean((res[,i] - res[,i - 1]) ^ 2))
  MSE <- mean((res[,i] - res[,i - 1]) ^ 2)

# variance of zSt
  VSE <- var((res[,i] - res[,i - 1]) ^ 2)

# 95% confidence using X-square distribution
# # lowwer limit
  ll.s <- qchisq(p = 0.025, df = N - 1, ncp = 0, lower.tail = TRUE, log.p = FALSE) * sqrt(VSE)
# # upper limit
  ul.s <- qchisq(p = 0.975, df = N - 1, ncp = 0, lower.tail = TRUE, log.p = FALSE) * sqrt(VSE)
  paste(round(MSE, 1), "ll:", round(ll.s, 1), "ul:", round(ul.s, 1))

############################################ SS (Sum of squares)
  SS <- sum((res[,i] - res[,i - 1]) ^ 2)

# fill report table
  report[i / 2,1:4] <- c(names(res)[i - 1], ME, RMSE, SS)
}

#d.stat <- read.csv("summary.calibdata.csv")
report$R2 <- 1 - (as.numeric(report[,4]) / d.stat[6,1:7])

#########Plot
library(semPlot)
#matrix to arrange nodes at graph
layout <-as.matrix(read.csv("matrix_semplot.csv")[,-1])
# plot sem using layout=layout

layout[]

semPaths(fit3,what = "std",whatLabels = "no", layout = layout,sizeLat = 4, cut =0.35,
         sizeInt2 = 2,sizeMan =3.5, sizeLat2 = 2, sizeInt = 1,sizeMan2 = 1.5,nCharNodes=0, font=3,
         edge.width = 2,esize=1.5, asize=1,intercepts = F, reorder = F,equalizeManifests =T, residuals = F,layoutSplit=F,
         structural = F, exoCov = F, exoVar=F,cardinal = F,style = "lisrel",#color = c("orange","blue"),
         manifests = c("tb.A", "sat.A", "evim", "evisd", "lstm", "lstsd", "dem","bt", "wdist", "mrvbf", "vdchn", 
                       "twi", "river", "thick.A", "slope", "maxc", "oc.A", "esp.B", "esp.A"))

semPlotModel(fit3)
semPaths(fit3)

str(fit3)

write.csv(report, "/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/report_cross-validation.csv")


semPaths




