############### #### ### ## # - ALTERNATIVE WAY TO DOWNLOAD MODIS DATA - # ## ### #### ###############
####                     #### For MOD13Q1 product EVI 16-days 250 m res. ####

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


rm(list=ls())
#setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/MODIS") # @RStudio desktop
setwd("~/big/USA/MODIS/") # @RStudio server
#write("/home/mangelini/big/temp/'", file=file.path(Sys.getenv('R_SESSION_TMPDIR'), '.Renviron')) # @RStudio
#library(raster)
library(rgdal)
#library(gdalUtils)
library(doParallel)

### from REVERB it is posible to get a csv file with a query result
# http://reverb.echo.nasa.gov/reverb/#utf8=%E2%9C%93&spatial_map=satellite&spatial_type=rectangle
# it can be directly read from ftp://ladsweb.nascom.nasa.gov/allData/5/MOD13Q1/ (Curl package)
# be aware of repeated files
files<-read.csv("MOD13Q1_reverb.csv") 

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
tile3 <- unique(granuleID$tile)[3]

# create directory to save files per tiles
# dir.create("output/")
# dir.create(path = paste("output/", unique(granuleID$tile)[1],"/", sep= ""),showWarnings = T)
# dir.create(path = paste("output/", unique(granuleID$tile)[2],"/", sep= ""),showWarnings = T)
# dir.create(path = paste("output/", unique(granuleID$tile)[3],"/", sep= ""),showWarnings = T)

### be careful, big data downloading 
# download data from tile j
tiles <- unique(granuleID$tile)
registerDoParallel(cores=10) # each core takes a i value
foreach(j = 1:3){
   n <- as.numeric(rownames(granuleID[granuleID$tile == unique(granuleID$tile)[j],]))
   # parallel processing 6 cores with doParallel package (it allows to open several downloads at once)
   # miss-connection could be solved with try() function
   foreach(i= n)  %dopar%{
      download.file(files$Online.Access.URLs[i],
                    destfile = paste("output/",granuleID$tile[i],"/", granuleID$tile[i],"_",
                                     granuleID$date[i], ".hdf", sep=""), quiet = TRUE,
                    mode = "wb", method = "wget", 
                    extra = "--load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies")
   }
}

# https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget
# mode="wb", method = "wget", extra = "--password=Icct5mily --user=angelini75")

# some issues with HDF4 images http://markmail.org/thread/thlchxcqf364wl5p
# extracting MODIS_Grid_16DAY_250m_500m_VI:250m 16 days EVI from HDF file
# %dopar% does not work here 
#foreach(i= 1)  %do%{

hdf <- paste0("output/",tiles[1],"/",list.files(paste0("~/big/USA/MODIS/output/",tiles[1]), ".hdf"))

hdfImage <- list()
foreach(i= seq_along(hdf))  %dopar%{
   hdfImage <- list()
   hdfImage[[i]] <- readGDAL(paste("HDF4_EOS:EOS_GRID:",hdf[i] , ":MODIS_Grid_16DAY_250m_500m_VI:250m 16 days EVI", sep = ""))
   # convert HDF to GeoTIFF
   writeGDAL(hdfImage[[i]], gsub(hdf[i], pattern = ".hdf",replacement= ".tif"),
             drivername = "GTiff", type = "Float32", mvFlag = NA, options=NULL, copy_drivername = "GTiff", setStatistics=FALSE) 
   # delete data to prevent large list
   hdfImage <- ""
   # delete HDF file to save space at disk
   unlink(hdf[i],force = T)
}

rm(hdf)
hdf <- paste0("output/",tiles[2],"/",list.files(paste0("~/big/USA/MODIS/output/",tiles[2]), ".hdf"))

hdfImage <- list()
foreach(i= seq_along(hdf))  %dopar%{
   hdfImage <- list()
   hdfImage[[i]] <- readGDAL(paste("HDF4_EOS:EOS_GRID:",hdf[i] , ":MODIS_Grid_16DAY_250m_500m_VI:250m 16 days EVI", sep = ""))
   # convert HDF to GeoTIFF
   writeGDAL(hdfImage[[i]], gsub(hdf[i], pattern = ".hdf",replacement= ".tif"),
             drivername = "GTiff", type = "Float32", mvFlag = NA, options=NULL, copy_drivername = "GTiff", setStatistics=FALSE) 
   # delete data to prevent large list
   hdfImage <- ""
   # delete HDF file to save space at disk
   unlink(hdf[i],force = T)
}

rm(hdf)
hdf <- paste0("output/",tiles[3],"/",list.files(paste0("~/big/USA/MODIS/output/",tiles[3]), ".hdf"))

hdfImage <- list()
foreach(i= seq_along(hdf))  %dopar%{
   hdfImage <- list()
   hdfImage[[i]] <- readGDAL(paste("HDF4_EOS:EOS_GRID:",hdf[i] , ":MODIS_Grid_16DAY_250m_500m_VI:250m 16 days EVI", sep = ""))
   # convert HDF to GeoTIFF
   writeGDAL(hdfImage[[i]], gsub(hdf[i], pattern = ".hdf",replacement= ".tif"),
             drivername = "GTiff", type = "Float32", mvFlag = NA, options=NULL, copy_drivername = "GTiff", setStatistics=FALSE) 
   # delete data to prevent large list
   hdfImage <- ""
   # delete HDF file to save space at disk
   unlink(hdf[i],force = T)
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
setwd("~/big/USA/MODIS")
startdir <- getwd()
library(rgdal)
library(gdalUtils)
library(raster)
library(doParallel)
library(maptools)
############ Create mosaic 
# subset extension in MODIS coordinate system (could be calculated from a shape file)
extent(readShapePoints("K_3_MODproj.shp"))

mask.ext <- c(-8995835.65,4052136.73, -7932725.14,4469948.27) # xmin,ymin , xmax,ymax @QGIS K_3MODproj.shp

# list of files to be processed
setwd("output")
tmp.lst <- list.files(include.dirs = T,full.names = T,pattern = ".tif$", recursive = T)
#raster projection of files (not needed yet)
MODIS <- projection(readGDAL(tmp.lst[1]))

# mosaiking is by pairs of tiles. The list is order by tiles, so after the last date from tile 1
# is the firs date of tile 2. (it is a bit risky step, could be safer <- improve!)
n <-length(tmp.lst)/3
#nombre de mosaicos
files.vrt <- NULL
for(i in 1:n){
  files.vrt<- append(x = files.vrt,gsub(".tif", ".vrt", gsub("h10v05/h10v05_A","",tmp.lst[i])))
}
# mosaiking
registerDoParallel(cores=10)
# log <- list() # optional
 foreach(i = 1:n) %dopar% {
   tmp.lst.i <-tmp.lst[c(i,i+346,i+346*2)]
   #write(files.vrt[1],file = files.vrt[i])
   gdalwarp(srcfile=tmp.lst.i, dstfile=files.vrt[i],te=mask.ext) # it mosaics and crops data with "te=" great!
#   log[[i]]<-  append(print(gdalinfo(files.vrt[i]))) # optional
 }
plot(raster("2008129.vrt"))
require(RCurl)
##read directories from ftp site (fixed MODIS periods). Function developed by Ype.
readdir<-function(fil){
  filenames <- try(getURL(fil,ftp.use.epsv = FALSE, dirlistonly = TRUE))
  filenames <- paste(fil, strsplit(filenames, "\r*\n")[[1]], sep = "")
  files <- filenames[grepl('\\.', substr(filenames,30,nchar(filenames)) )|
                       grepl('readme', substr(filenames,30,nchar(filenames)) )|
                       grepl('README', substr(filenames,30,nchar(filenames)) )]
  dirs <- setdiff(filenames, files)
  return(dirs)
}
# get MODIS periods
period <- gsub("ftp://ladsweb.nascom.nasa.gov/allData/5/MOD13Q1/2005/","",
                readdir("ftp://ladsweb.nascom.nasa.gov/allData/5/MOD13Q1/2005/"))
period
# in case I need projection from my study area I read this (tiff)image
#t_pj <- projection(readGDAL("./h12v12/target_pj"))

# crate list of images from the same period
pb = txtProgressBar(min = 0, max = length(period), initial = 0, style = 3) # see progress 1/2
registerDoParallel(cores=10) 
foreach(k = seq_along(period)) %dopar%{
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
