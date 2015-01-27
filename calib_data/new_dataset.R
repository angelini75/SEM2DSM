setwd("/media/L0135974_DATA/UserData/BaseARG/2_Calibration")

D <-read.csv("calib.data-0.1.1.csv")

names(D)
d.ph <- cbind(D[,c(1:4,94,95,98,7:18,20:26,30:33,51,55,61,77,81,87,88,90,91)])
as.numeric(d.ph$a_ph_kcl == 73.2) 
as.numeric(d.ph$a_ph_kcl == 0.3) 
d.ph[1482,12]<- 7.32
d.ph[1516,12]<- NA
boxplot(d.ph$a_ph_h2o-d.ph$a_ph_kcl)
