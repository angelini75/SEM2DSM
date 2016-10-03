############### #### ### ## # - ALTERNATIVE WAY TO DOWNLOAD MODIS DATA - # ## ### #### ###############
####                     #### For MCD43A4 product reflectance bands 500 m res. ####

# Purpose        : Download 15 years tiles MODIS product (MOD13Q1); mosaic them; substract AOI; calculate mean and sd per period; 
# Maintainer     : Marcos Angelini  (marcos.angelini@wur.nl); 
# Contributions  : Tomislav Hengl; Ype van der Velde; Wei Shangguan
# Status         : pre-pre-alpha
# Note           : 
# sessionInfo(@RStudio desktop)  lenovo ThinkPad T430 (4 cores)
# R version 3.0.2 (2013-09-25)
# Platform: x86_64-pc-linux-gnu (64-bit)
# sessionInfo(@RStudio Server) http://h2335862.stratoserver.net:8787/ (12 cores)
# R version 3.1.1 (2014-07-10)
# Platform: x86_64-unknown-linux-gnu (64-bit)
# fileh12v12 <- list.files(path = "h12v12/",pattern = ".hdf$",)
# fileh13v12 <- list.files(path = "h13v12/",pattern = ".hdf$",)
# file1 <- cbind(tile = "h12", file = fileh12v12, size = file.size(paste0(path = "h12v12/",fileh12v12))/1000000)
# file2 <- cbind(tile = "h13",file = fileh12v12, size = file.size(paste0(path = "h12v12/",fileh12v12))/1000000)
# file3 <- as.data.frame(rbind(file1,file2))
# file3$size <- as.character(file3$size)  
# file3$size <- as.double(file3$size)  
# miss <- as.character(file3$file[file3$size<20])
# files<-read.csv("save_results_csv.csv") 
# name(files)
# files <- files[files$Producer.Granule.ID == miss,]
# 
# rm(list=ls())
# setwd("/home/marcos/MCD43A4") # @RStudio desktop
# #setwd("~/modis") # @RStudio server
# #library(raster)
# library(rgdal)
# #library(gdalUtils)
# library(doParallel)
# name <- function(x) { as.data.frame(names(x))} 
# 
# ### from REVERB it is posible to get a csv file with a query result
# # http://reverb.echo.nasa.gov/reverb/#utf8=%E2%9C%93&spatial_map=satellite&spatial_type=rectangle
# # it can be directly read from ftp://ladsweb.nascom.nasa.gov/allData/5/MOD13Q1/ (Curl package)
# # be aware of repeated files
# 

setwd(dir = "sftp://mangelini@85.214.221.253/home/mangelini/bigdir")

f <-list.files("output/h12v12/", pattern = ".hdf$")
f <- gsub(pattern = "h12v12_",replacement = "",x = f)
mask.ext <- raster(xmn=-5672223,ymn=-3899255,xmx= -5384129,ymx= -3693436,nrow=100,ncol=100)
mask.ext[] <- runif(100*100)
pj <- projection(hdfImage[[1]])
pb <- txtProgressBar(min = 0, max = length(f), initial = 0, style = 3)

for(i in 1:length(f)){
hdfImage <- list() 
hdfImage[[1]] <- readGDAL(paste0("HDF4_EOS:EOS_GRID:","output/h12v12/","h12v12_",f[i], ":MOD_Grid_BRDF:Nadir_Reflectance_Band3"))
hdfImage[[2]] <- readGDAL(paste0("HDF4_EOS:EOS_GRID:","output/h13v12/","h13v12_",f[i], ":MOD_Grid_BRDF:Nadir_Reflectance_Band3"))
n <- merge(raster(hdfImage[[1]]),raster(hdfImage[[2]]))
m <- crop(x = n,y = mask.ext)
raster::writeRaster(x = m, filename = paste0("/home/marcos/MCD43A4/B3/",f[i],".B3.tif"), overwrite=TRUE)
setTxtProgressBar(pb,i)
}
# create a random raster over the space:        


# plot it with the boundaries we want to clip against:
plot(n)
plot(mask.ext,add=TRUE)

# now use the mask function
rr <- mask(mask = mask.ext,x =  n)
plot(raster::)
# convert HDF to GeoTIFF
  writeGDAL(hdfImage[[1]], paste("output/TIFF", f[1],"B3", ".tif", sep=""),
            drivername = "GTiff", type = "Float32", mvFlag = NA, options=NULL, copy_drivername = "GTiff", setStatistics=FALSE) 
  
  # delete data to prevent large list
  hdfImage <- ""
  # delete HDF file to save space at disk
  unlink(paste("output/",granuleID$tile[i],"/",
               granuleID$tile[i],"_",granuleID$date[i], ".hdf", sep=""),force = T)
}
}

#### END DOWNLOADING!

## control 
# filelist <- list.files(path = paste("output/",unique(granuleID$tile)[1],"/",  sep=""),pattern = "tif")
# sum(as.numeric(granuleID$tile == "h13v12")) > length(filelist) 
# target <- paste(unique(granuleID$tile)[1], "_",granuleID$date[granuleID$tile == unique(granuleID$tile)[1]], ".tif", sep="")
# setdiff(target, filelist)
# #   b_tmp[b_tmp%in%a_tmp]   
# length(target %in% filelist)


############### #### ### ## # - IMAGE PROCESSING - # ## ### #### ###############
rm(list=ls())
setwd("/media/L0135974_DATA/UserData/BaseARG/COVARIATES/MODIS")
startdir <- getwd()
library(rgdal)
library(gdalUtils)
library(raster)
library(doParallel)
############ Create mosaic 
# subset extension in MODIS coordinate system (could be calculated from a shape file)
mask.ext <- c(-5672223,-3899255, -5384129, -3693436)

# list of files to be processed
setwd("output")
tmp.lst <- list.files(include.dirs = T,full.names = T,pattern = ".tif$", recursive = T)
#raster projection of files (not needed yet)
pj <- projection(readGDAL(tmp.lst[1]))

# mosaiking is by pairs of tiles. The list is order by tiles, so after the last date from tile 1
# is the firs date of tile 2. (it is a bit risky step, could be safer <- improve!)
n <-length(tmp.lst)/2
#nombre de mosaicos
files.vrt <- NULL
for(i in 1:n){
  files.vrt<- append(x = files.vrt,gsub(".tif", ".vrt", gsub("h12v12/h12v12_A","",tmp.lst[i])))
}
# mosaiking
registerDoParallel(cores=10)
# log <- list() # optional
 foreach(i = 1:n) %dopar% {
   tmp.lst.i <-tmp.lst[c(i,i+length(tmp.lst)/2)]
   #write(files.vrt[1],file = files.vrt[i])
   gdalwarp(srcfile=tmp.lst.i, t_srs=pj, dstfile=files.vrt[i],te=mask.ext) # it mosaics and crops data with "te=" great!
#   log[[i]]<-  append(print(gdalinfo(files.vrt[i]))) # optional
 }

require(RCurl)
##read directories from ftp site (fixed MODIS periods). Function developed by Ype.
readdir<-function(fil){
  filenames <- try(getURL(fil,ftp.use.epsv = FALSE, dirlistonly = TRUE))
  filenames <- paste(fil, strsplit(filenames, "\r*\n")[[1]], sep = "")
  files <- filenames[grepl('\\.', substr(filenames,30,nchar(filenames)) )|grepl('readme', substr(filenames,30,nchar(filenames)) )|grepl('README', substr(filenames,30,nchar(filenames)) )]
  dirs <- setdiff(filenames, files)
  return(dirs)
}
# get MODIS periods
period <- gsub("ftp://ladsweb.nascom.nasa.gov/allData/5/MOD13Q1/2005/","",
               readdir("ftp://ladsweb.nascom.nasa.gov/allData/5/MOD13Q1/2005/"))

# in case I need projection from my study area I read this (tiff)image
#t_pj <- projection(readGDAL("./h12v12/target_pj"))

# crate list of images from the same period
pb = txtProgressBar(min = 0, max = length(period), initial = 0, style = 3) # see progress 1/2
for(k in 1:length(period)){
  list.per <- paste("./", list.files(pattern= paste(period[k],".vrt",sep="")),sep="")
  a <-list()
# stack images in a list (not needed)
#   for (j in 1:length(list.per)) {
#     a[[j]] <- readGDAL(list.per[j])@data$band1/10000 # redundant 
#     a[[j]][a[[j]]<0] <-NA # clean negative values
#   }
  setTxtProgressBar(pb,k) # see progress 2/2
# stack images from the same period (raster package)
# I found this script here http://markmail.org/thread/lsjnczi5fppddrpr
  s <- stack(list.per)
  s <- s/10000 # to reduce no significan digits
  print(paste("###############   number of layer period ",period[k],": ",nlayers(x = s),sep="")) # just FYI
# to calculate the mean of stacked rasters (very simple!)
  m<-mean(s,na.rm=T) 
  #plot(m) # if you want
# To calculate sd (because sd(s) got error)
  std<-calc(s, fun = sd)
  #plot(std) # if you want

    # Optionally, you can reproject and save the raster with gdalwarp() 
    # m.sd <- as(m,"SpatialPixelsDataFrame")
    # gdalwarp(srcfile=m.sd, s_srs=pj,t_srs=t_pj, dstfile= paste(period[k],"m.tif"),tr=c(250,250),r="bilinear")

# save rstarers (mean and sd) as tiff file. Each per MODIS period. 
# http://r-sig-geo.2731867.n2.nabble.com/Problems-with-writeGDAL-gridded-function-td5385974.html
# writeGDAL can't deal with raster, for this reason writeRaster() exists. 
  writeRaster(x = m,filename =  paste("./",period[k],"_mean.tif",sep=""), overwrite=T)
  writeRaster(x = std,filename =  paste("./",period[k],"_std.tif",sep=""), overwrite=T)
  
  rm(s)
  rm(m)
  rm(std)
}

#### End processing (so far)
