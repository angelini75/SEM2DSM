# This script is to produce the results shown in Wageningen Soil Conference presentation
# Before to run it, it is necessary to run First_model_v0.2.R til prediction section


###### OC prediction and comparison between SEM and linear regression
# sem oc
library(utils)
pb = txtProgressBar(min = 0, max = length(pred[,1]), initial = 0, style = 3)
pred <- pred[,15:29]
for(i in 1:length(pred[,1])) {
  p=matrix(c(pred$evim[i], pred$evisd[i],
             pred$lstm[i], pred$lstsd[i], pred$dem[i], 
             pred$wdist[i], pred$mrvbf[i], pred$vdchn[i],
             pred$twi[i], pred$river[i], pred$slope[i],
             pred$maxc[i]),nrow=12,ncol=1)
  n=c(0,0,0,0,0,0,0)
  n=(solve(I-B))%*%((A%*%p)) # key equation
  #   pred$tb.Ar[i] <-  n[1]
  #   pred$sat.Ar[i] <-  n[2]
  #   pred$btr[i] <-  n[3]
  pred$oc.Ar[i] <-  n[4]
  #   pred$thick.Ar[i] <-  n[5]
  #   pred$esp.Br[i] <-  n[6]
  #   pred$esp.Ar[i] <-  n[7]
  setTxtProgressBar(pb,i)
}

# lm oc
oc.fit<- lm(formula = oc.A ~ lstm +  lstsd + evim + evisd + dem + wdist + mrvbf + vdchn + twi, data = d)
summary(oc.fit)
pred$oc.lm <-as.vector(predict(oc.fit,pred))


as.data.frame(names(pred))
write.csv(pred, "pred.bk.oc.csv")
pred <- pred[,14:17]
N <- N[c(4),]

pred <- read.csv("pred.bk.oc.csv")
pred <- pred[,-1]
pred$oc.Ar<- pred$oc.Ar*N[6,3] + N[6,2]
pred$oc.lm<- pred$oc.lm*N[6,3] + N[6,2]
as.data.frame(names(pred))
pred <- pred[,14:17]

library(sp)
library(raster)
pred.sp <- pred
coordinates(pred.sp) <- ~X+Y
# spplot(pred.sp)

y <- raster("mask_231m_posgar.tif")
proj4string(y) <- posgar98
proj4string(pred.sp) <- posgar98
#pred.sp <- spTransform(pred.sp, modis)
r <- rasterize(x = pred.sp,y = y,background= NA)
plot(r)
writeRaster(x = r,filename ="rusults.tif", overwrite=T,bylayer=TRUE,suffix=r@data@names)

summary(oc.lm)

