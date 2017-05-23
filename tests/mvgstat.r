# script taken from a cokriging study
# script below defines variograms and cross-variograms for 6 variables.
# see also http://www.sciencedirect.com/science/article/pii/S2352009416300219
setwd("/home/marcos/Documents/SEM2DSM1/tests")
library(sp)
library(gstat)

# read observations:
reslogOC = read.csv("logOCres.csv", header = TRUE, sep = ",", na.string = "-999")
reslogTOTN = read.csv("logTOTNres.csv", header = TRUE, sep = ",", na.string = "-999")
reslogOT = cbind(reslogOC,reslogTOTN$logTOTNAres,
                 reslogTOTN$logTOTNBres,reslogTOTN$logTOTNCres)

# add noise to coordinates to avoid duplicates 
reslogOT$X = reslogOT$X + rnorm(dim(reslogOT)[1], mean = 0, sd = 10)
reslogOT$Y = reslogOT$Y + rnorm(dim(reslogOT)[1], mean = 0, sd = 10)

# rename variables
names(reslogOT)[6:11] = c("logOC_A","logOC_B","logOC_C","logON_A","logON_B","logON_C")
names(reslogOT)

# copy data and create a spatial data set
spreslogOT = reslogOT
coordinates(spreslogOT) = ~X+Y
proj4string(spreslogOT) = CRS("++init=epsg:3035")

# variography requires that all variables are measured at all points
vgreslogOT = subset(reslogOT,!is.na(reslogOT[6]) & !is.na(reslogOT[7]) 
                    &!is.na(reslogOT[8]) & !is.na(reslogOT[9]) & !is.na(reslogOT[10]) 
                    &!is.na(reslogOT[11]))
coordinates(vgreslogOT) = ~X+Y
proj4string(vgreslogOT) = CRS("++init=epsg:3035")

# define parameters of Linear Model of Coregionalization
Nv = 6
range = 1800000
nugget = c(0.57, 0.55, 0.7, 0.48, 0.68, 0.85)
psill = 1 - nugget
rho = array(1, dim = c(Nv,Nv))
rho[2,1]     = c(0.42)
rho[3,1:2]   = c(0.4,0.62)
rho[4,1:3]   = c(0.82,0.34,0.35)
rho[5,1:4]   = c(0.41,0.8,0.56,0.42)
rho[6,1:5]   = c(0.34,0.45,0.75,0.4,0.65)
for (i in 1:(Nv-1)) {
  for (j in (i+1):Nv) {
    rho[i,j] = rho[j,i] }  }
# eigen(rho,only.values=TRUE)

# define multivariate gstat object 
glogOT = gstat(id = c("C_A"), formula = logOC_A~1, data = vgreslogOT,
               beta = 0, model = vgm(psill[1],"Sph",range,nugget[1]))
glogOT = gstat(glogOT, id = "C_B", formula = logOC_B~1, data = vgreslogOT,
               beta = 0, model = vgm(psill[2],"Sph",range,nugget[2]))
glogOT = gstat(glogOT, id = "C_C", formula = logOC_C~1, data = vgreslogOT,
               beta = 0, model = vgm(psill[3],"Sph",range,nugget[3]))
glogOT = gstat(glogOT, id = "N_A", formula = logON_A~1, data = vgreslogOT,
               beta = 0, model = vgm(psill[4],"Sph",range,nugget[4]))
glogOT = gstat(glogOT, id = "N_B", formula = logON_B~1, data = vgreslogOT,
               beta = 0, model = vgm(psill[5],"Sph",range,nugget[5]))
glogOT = gstat(glogOT, id = "N_C", formula = logON_C~1, data = vgreslogOT,
               beta = 0, model = vgm(psill[6],"Sph",range,nugget[6]))

names = c("C_A","C_B","C_C","N_A","N_B","N_C")
for (i in 1:(Nv-1)) {
  for (j in (i+1):Nv) {
    glogOT = gstat(glogOT, id=c(names[i],names[j]), 
                   model=vgm(rho[i,j]*sqrt(psill[i]*psill[j]),
                             "Sph",range,rho[i,j]*sqrt(nugget[i]*nugget[j]))) } }
windows(width = 12, height = 12)
plot(variogram(glogOT), model = glogOT$model, main = "")
savePlot(filename = "Figure 9 vgmlogOCTOTN", type = "png")
