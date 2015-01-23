setwd("/media/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")

D <- read.csv("calib.data-1.0.csv")

# CEC from OC (Nyle and Weil, 2007: chapter 8)
D$CEC_OM <- (D$a_ph_h2o*50/1.5-33.33)*(D$a_mo_c*1.72/100)
# CEC from clay + silt
D$CEC_cl.sl <- D$a_t - D$CEC_OM 
# CEC from silt considering CEC of silt = 15.5 cmol/kg (range 8 to 23 from Morras, 1995)
D$CEC_sl <- 15.5*(D$a_limo_2_2/100)
# CEC from clay
D$CEC_cl <- D$a_t - D$CEC_OM - D$CEC_sl
# CEC from pure clay 
D$CEC_pure.cl <- D$CEC_cl/(D$a_arcilla/100)


par(mfrow = c(2, 2))
a <- c("A", "Bt", "BC", "C")
for(i in a){
  hist(D$CEC_pure.cl[D$hor == i],breaks = 20, xlab = paste("Clay CEC ",i," horizon",sep=""))
}

# spatial distribution
library(sp)
d <- cbind(D[D$hor == "BC",1:2],CEC_cl=D$CEC_pure.cl[D$hor == "BC"])
d <- na.omit(d)
d$CEC_cl <- round(d$CEC_cl,0)
coordinates(d) <- ~X+Y
spplot(d,zcol ="CEC_cl",edge.col="black", colorkey=T,col.regions= rainbow(1000) )

# Morrás, Héctor JM. "Mineralogy and cation exchange capacity of the fine silt fraction in 
# two soils from the southern Chaco Region (Argentina)." Geoderma 64, no. 3 (1995): 281-295.

# Brady, Nyle C., and Ray R. Weil. "The nature and properties of soils." No. Ed. 14.
# Prentice-Hall Inc., 2007.

##### history()
# plot(D$a_t,D$a_mo_c)
# plot(D$a_mo_c, D$a_t)
# plot(D$a_mo_c, D$a_t)[D$horizon == "A"]
# plot(D$a_mo_c, D$a_t)[D$hor == "A"]
# plot(D$a_mo_c[D$hor == "A"], D$a_t)
# plot(D$a_mo_c[D$hor == "A"], D$a_t[D$hor == "A"])
# plot(D$a_arcilla[D$hor == "A"], D$a_t[D$hor == "A"])
# plot((D$a_arcilla[D$hor == "A"])/(D$a_mo_c[D$hor == "A"]), D$a_t[D$hor == "A"])
# plot((D$a_arcilla[D$hor == "A"])*(D$a_mo_c[D$hor == "A"]), D$a_t[D$hor == "A"])
# plot(D$a_arcilla[D$hor == "A"], D$a_t[D$hor == "A"])
# plot(D$a_arcilla[D$hor == "Bt"], D$a_t[D$hor == "Bt"])
# plot(D$a_arcilla[D$hor == "A"], D$a_mo_c[D$hor == "A"])
# plot(D$a_arcilla[D$hor == "E"], D$a_t[D$hor == "E"])
# plot(D$a_mo_c[D$hor == "A"], D$a_t[D$hor == "A"])
# 50/1.5*5.5
# (50/1.5*5.5)+100
# (1.5/50*5.5)+100
# (50/1.5*5.5)+100
# (50/1.5*5.5)+10
# 50/1.5
# 50/1.5*4
# 50/1.5*4-33
# D <- read.csv("calib.data-1.0.csv")
# (D$a_ph_h2o*50/1.5-33.33)*(D$a_mo_c*1.72/100)
# boxplot(D$a_ph_h2o*50/1.5-33.33)*(D$a_mo_c*1.72/100)
# (D$a_mo_c*1.72/100)
# boxplot(D$a_mo_c*1.72/100)
# (D$a_ph_h2o*50/1.5-33.33)
# (7*50/1.5-33.33)
# 200*.04
# (D$a_ph_h2o*50/1.5-33.33)*(D$a_mo_c*1.72/100)
# boxplot(D$a_ph_h2o*50/1.5-33.33)*(D$a_mo_c*1.72/100)
# (D$a_ph_h2o*50/1.5-33.33)*(D$a_mo_c*1.72/100)
# range(D$a_ph_h2o*50/1.5-33.33)*(D$a_mo_c*1.72/100)
# mean(D$a_ph_h2o*50/1.5-33.33)*(D$a_mo_c*1.72/100)
# boxplot((D$a_ph_h2o*50/1.5-33.33)*(D$a_mo_c*1.72/100))
# D$CEC_OM <- (D$a_ph_h2o*50/1.5-33.33)*(D$a_mo_c*1.72/100)
# D$CEC_Clay <- D$a_t - D$CEC_OM
# summary(D$CEC_Clay)
# boxplot(D$CEC_Clay)
# D$CEC_Clay*(D$a_arcilla/100)
# hist(D$CEC_Clay*(D$a_arcilla/100))
# (D$a_arcilla/100)
# hist(D$CEC_Clay/(D$a_arcilla/100))
# hist(D$CEC_Clay/(D$a_arcilla/100),freq = T)
# hist(D$CEC_Clay/(D$a_arcilla/100),breaks = 10)
# hist(D$CEC_Clay/(D$a_arcilla/100),breaks = 2)
# hist(D$CEC_Clay/(D$a_arcilla/100),breaks = 20)
# hist(D$CEC_Clay/(D$a_arcilla/100),breaks = 50)
# CEC_pure.cl <- D$CEC_Clay/(D$a_arcilla/100)
# D$CEC_pure.cl <- D$CEC_Clay/(D$a_arcilla/100)
# plot(D$CEC_pure.cl, D$a_limo_2_5)
