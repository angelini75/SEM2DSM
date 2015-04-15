############### #### ### ## # - ALTERNATIVE WAY TO DOWNLOAD MODIS DATA - # ## ### #### ###############
####                     #### For MOD11A2  ####
#    Terra Land Surface Temperature & Emissivity | Tile 1000m | 8 day 
# Purpose        : Download 15 years tiles MODIS product (MOD11A2); mosaic them; substract AOI; calculate mean and sd per period; 
# Maintainer     : Marcos Angelini  (marcos.angelini@wur.nl); 
# Contributions  : 
# Status         : on going
# Note           : 
# sessionInfo(@RStudio desktop)  lenovo ThinkPad T430 (4 cores)
# R version 3.0.2 (2013-09-25)
# Platform: x86_64-pc-linux-gnu (64-bit)
# sessionInfo(@RStudio Server) http://h2335862.stratoserver.net:8787/ (12 cores)
# R version 3.1.1 (2014-07-10)
# Platform: x86_64-unknown-linux-gnu (64-bit)


rm(list=ls())
#setwd("/media/L0135974_DATA/UserData/BaseARG/COVARIATES/MODIS") # @RStudio desktop
setwd("~/big") # @RStudio server
#library(raster)
library(rgdal)
#library(gdalUtils)
library(doParallel)

### from REVERB it is posible to get a csv file with a query result
# http://reverb.echo.nasa.gov/reverb/#utf8=%E2%9C%93&spatial_map=satellite&spatial_type=rectangle
# it can be directly read from ftp://ladsweb.nascom.nasa.gov/allData/5/MOD13Q1/ (Curl package)
# be aware of repeated files
files<-read.csv("save_results_csv_MOD11A2.csv") 

# to be aware
paste(length(files$Granule.UR)," files here!", sep="")
paste(sum(as.numeric(gsub(",",".",files$Size)))/1000, "Gb!")

# URL of each HDF4 file
files$id <- 1:length(files[,1])
files$Online.Access.URLs <- as.character(files$Online.Access.URLs)
granule <- files[,c(2,10)]
write.table(granule, "granule.txt",row.names = F,col.names = F,quote = F, sep=".")

# espesifications of each file
granuleID <- read.csv("granule.txt", sep=".", header=F)
names(granuleID) <- c("product", "date","tile","ver","proj","format","id")

## Tiles
granuleID$date <- as.character(granuleID$date)
granuleID$tile <- as.character(granuleID$tile)
tile1 <- unique(granuleID$tile)[1]
tile2 <- unique(granuleID$tile)[2]

# create directory to save files per tiles
dir.create("output/")
dir.create(path = paste("output/", unique(granuleID$tile)[1],"/", sep= ""),showWarnings = T)
dir.create(path = paste("output/", unique(granuleID$tile)[2],"/", sep= ""),showWarnings = T)

### be careful, big data downloading 
# download data from tile j
registerDoParallel(cores=6) # each core takes a i value
for (j in 1:2){
  n <- as.numeric(rownames(granuleID[granuleID$tile == unique(granuleID$tile)[2],]))
  # parallel processing 6 cores with doParallel package (it allows to open several downloads at once)
  # miss-connection could be solved with try() function
  hdfImage <- list()
  
  foreach(i = n)  %dopar%{
    download.file(files$Online.Access.URLs[i],destfile = paste("output/",granuleID$tile[i],"/",
                                                               granuleID$tile[i],"_",granuleID$date[i], ".hdf", sep=""), mode="wb")
  }
}

f <-list.files("output/h12v12/", pattern = ".hdf")
  
    # some issues with HDF4 images http://markmail.org/thread/thlchxcqf364wl5p
    # extracting MODIS_Grid_16DAY_250m_500m_VI:250m 16 days EVI from HDF file
    # %dopar% does not work here 
foreach(i = 1:length(f))  %do%{
      hdfImage <- list()
      hdfImage[[i]] <- readGDAL(paste("HDF4_EOS:EOS_GRID:", paste("output/h12v12/", f[i], sep=""), 
                                      ":MODIS_Grid_8Day_1km_LST:LST_Day_1km", sep = ""))
    
    # convert HDF to GeoTIFF
    writeGDAL(hdfImage[[i]], paste("output/h12v12/", gsub(".hdf",".tif",f)[i], sep=""),
              drivername = "GTiff", type = "Float32", mvFlag = NA, options=NULL, copy_drivername = "GTiff", setStatistics=FALSE) 
}
# delete hdf files at h12v12 directory
for(i in 1:length(f)){
unlink(paste("output/h12v12/", f[i], sep=""),force = T)
}
# list hdf files
f <-list.files("output/h13v12/", pattern = ".hdf")
hdfImage <- list()
foreach(i = 1:length(f))  %do%{
  hdfImage <- list()
  hdfImage[[i]] <- readGDAL(paste("HDF4_EOS:EOS_GRID:", paste("output/h13v12/", f[i], sep=""), 
                                  ":MODIS_Grid_8Day_1km_LST:LST_Day_1km", sep = ""))

  writeGDAL(hdfImage[[i]], paste("output/h13v12/", gsub(".hdf",".tif",f)[i], sep=""),
            drivername = "GTiff", type = "Float32", mvFlag = NA, options=NULL, copy_drivername = "GTiff", setStatistics=FALSE) 
  hdfImage <- ""
}
# delete hdf files at h13v12 directory
for(i in 1:length(f)){
  unlink(paste("output/h13v12/", f[i], sep=""),force = T)
}

#  }
gc(verbose = getOption("verbose"), reset = T)
#}

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
setwd("~/big") # @RStudio server
# setwd("/media/L0135974_DATA/UserData/BaseARG/COVARIATES/MODIS")
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

# mosaiking is made by pairs of tiles. The list is order by tiles, so after the last date from tile 1
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
period <- gsub("ftp://ladsweb.nascom.nasa.gov/allData/5/MOD11A2/2005/","",
               readdir("ftp://ladsweb.nascom.nasa.gov/allData/5/MOD11A2/2005/"))

# in case I need projection from my study area I read this (tiff)image
#t_pj <- projection(readGDAL("./h12v12/target_pj"))

# crate list of images from the same period
# pb <- txtProgressBar(min = 0, max = length(period), initial = 0, style = 3) # see progress 1/2
# for(k in 1:length(period)){
  # list.per <- paste("./", list.files(pattern= paste(period[k],".vrt",sep="")),sep="")
  list.tot <- paste("./", list.files(pattern= ".vrt"),sep="")
  a <-list()
  # stack images in a list (not needed)
  #   for (j in 1:length(list.per)) {
  #     a[[j]] <- readGDAL(list.per[j])@data$band1/10000 # redundant 
  #     a[[j]][a[[j]]<0] <-NA # clean negative values
  #   }
  # setTxtProgressBar(pb,k) # see progress 2/2
  # stack images from the same period (raster package)
  # I found this script here http://markmail.org/thread/lsjnczi5fppddrpr
  s <- stack(list.tot)
  #s <- s/10000 # to reduce no significan digits
  #print(paste("###############   number of layer period ",period[k],": ",nlayers(x = s),sep="")) # just FYI
  # to calculate the mean of stacked rasters (very simple!)
  m<-mean(s,na.rm=T) # times 0.02 to scale to kelvin (https://lpdaac.usgs.gov/products/modis_products_table/myd11a2)
  plot(m) # if you want
  # To calculate sd (because sd(s) got error)
  std<-calc(s, fun = sd, na.rm=T)
  plot(std) # if you want

  # Optionally, you can reproject and save the raster with gdalwarp() 
  # m.sd <- as(m,"SpatialPixelsDataFrame")
  # gdalwarp(srcfile=m.sd, s_srs=pj,t_srs=t_pj, dstfile= paste(period[k],"m.tif"),tr=c(250,250),r="bilinear")
  
  # save rstarers (mean and sd) as tiff file. 
  # http://r-sig-geo.2731867.n2.nabble.com/Problems-with-writeGDAL-gridded-function-td5385974.html
  # writeGDAL can't deal with raster, for this reason writeRaster() exists. 
#   writeRaster(x = m,filename =  paste("./",period[k],"_mean.tif",sep=""), overwrite=T)
#   writeRaster(x = std,filename =  paste("./",period[k],"_std.tif",sep=""), overwrite=T)
writeRaster(x = m,filename =  "mod11a2_tot_mean.tif", overwrite=T)
writeRaster(x = std,filename =  "mod11a2_tot_std.tif", overwrite=T)
#   rm(s)
#   rm(m)
#   rm(std)
#}

#### End processing (so far)

system('kill -15 19107')