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
rm(list=ls())
#setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/MODIS") # @RStudio desktop
setwd("~/big/USA/MODIS/")
library(doParallel)
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

readfiles<-function(fil){
  filenames <- try(getURL(fil,ftp.use.epsv = FALSE, ftplistonly = TRUE))
  filenames <- paste(fil, strsplit(filenames, "\r*\n")[[1]], sep = "")
  files <- filenames[grepl('//.', substr(filenames,30,nchar(filenames)) )|
                       grepl('readme', substr(filenames,30,nchar(filenames)) )|
                       grepl('README', substr(filenames,30,nchar(filenames)) )]
  filenames
}

#MODIS URL
URL <- "ftp://ladsweb.nascom.nasa.gov/allData/5"
#define MODIS product
MODISP <- "MCD43A4"
# driver and bands of HDF file for GDAL
driver <- "HDF4_EOS:EOS_GRID:"
band1 <- ":MOD_Grid_BRDF:Nadir_Reflectance_Band2"
band2 <- ":MOD_Grid_BRDF:Nadir_Reflectance_Band5"
#define tiles
tiles <- c("h09v05", "h10v04", "h10v05")
#define years
yrs <- as.character(2000:2015)
#get periods
period <- gsub("ftp://ladsweb.nascom.nasa.gov/allData/5/MCD43A4/2005/","",
               readdir("ftp://ladsweb.nascom.nasa.gov/allData/5/MCD43A4/2005/"))
# get file urls.
u <- NULL
for(i in seq_along(yrs)){
  for(j in seq_along(period)){
    u <-append(u,(paste0(URL,"/",MODISP,"/", yrs[i],"/",period[j], "/")))
  }
}

u <- u[-1:-6]

# TILE #1
n <- NULL 
for(i in seq_along(u)){
  z <- readfiles(u[i])
  n[i] <- z[grep(tiles[1], z)][1]
  ifelse(test = is.na(n[i]),
         n[i] <- z[grep(tiles[1], z)][1],
         n[i])
}
u1 <- u[which(is.na(n))]

# download HDF from urls
registerDoParallel(cores=10)
foreach(i = seq_along(n)) %dopar%{
  download.file(n[i], 
                destfile = paste("output/",tiles[1],"/",
                                 substr(x = n[i],start = 66,stop = 73),".hdf", sep=""), 
                quiet = TRUE, mode = "wb", method = "wget")
  }

# TILE #2
n <- NULL 
for(i in seq_along(u)){
  z <- readfiles(u[i])
  n[i] <- z[grep(tiles[2], z)][1]
  ifelse(test = is.na(n[i]),
         n[i] <- z[grep(tiles[2], z)][1],
         n[i])
}
u2 <- u[which(is.na(n))]

# download HDF from urls
registerDoParallel(cores=10)
foreach(i = seq_along(n)) %dopar%{
  download.file(n[i], 
                destfile = paste("output/",tiles[2],"/",
                                 substr(x = n[i],start = 66,stop = 73),".hdf", sep=""), 
                quiet = TRUE, mode = "wb", method = "wget")
}
u2 <- c("ftp://ladsweb.nascom.nasa.gov/allData/5/MCD43A4/2007/145/",
        "ftp://ladsweb.nascom.nasa.gov/allData/5/MCD43A4/2010/113/")


# TILE #3
n <- NULL 
for(i in seq_along(u)){
  z <- readfiles(u[i])
  n[i] <- z[grep(tiles[3], z)][1]
  ifelse(test = is.na(n[i]),
         n[i] <- z[grep(tiles[3], z)][1],
         n[i])
}
u3 <- u[which(is.na(n))]

# download HDF from urls
registerDoParallel(cores=10)
foreach(i = seq_along(n)) %dopar%{
  download.file(n[i], 
                destfile = paste("output/",tiles[3],"/",
                                 substr(x = n[i],start = 66,stop = 73),".hdf", sep=""), 
                quiet = TRUE, mode = "wb", method = "wget")
}

# rescuing missing files
granules <- read.csv(file = "MCD43A4_2002-2010-2007-2001.csv")
miss <- c("MCD43A4.A2002273.h09v05", "MCD43A4.A2007009.h09v05", "MCD43A4.A2007217.h09v05",
          "MCD43A4.A2007281.h09v05", "MCD43A4.A2007337.h09v05", "MCD43A4.A2010057.h09v05",
          "MCD43A4.A2007145.h10v04", "MCD43A4.A2010113.h10v04", "MCD43A4.A2001353.h10v05",
          "MCD43A4.A2001361.h10v05", "MCD43A4.A2007273.h10v05", "MCD43A4.A2007321.h10v05")

granules <- read.csv(file = "MCD43A4.csv")
miss2 <- c("MCD43A4.A2002073.h10v05", "MCD43A4.A2007305.h10v05", "MCD43A4.A2002033.h10v05",
           "MCD43A4.A2008049.h10v05", "MCD43A4.A2002137.h10v05", "MCD43A4.A2007241.h10v05",
           "MCD43A4.A2002025.h10v05", "MCD43A4.A2002065.h10v05", "MCD43A4.A2007353.h10v05",
           "MCD43A4.A2002113.h10v05", "MCD43A4.A2002129.h10v05", "MCD43A4.A2007249.h10v05",
           "MCD43A4.A2008009.h10v05", "MCD43A4.A2002217.h10v05", "MCD43A4.A2002177.h10v05",
           "MCD43A4.A2008001.h10v05", "MCD43A4.A2002105.h10v05", "MCD43A4.A2002185.h10v05",
           "MCD43A4.A2002065.h10v04", "MCD43A4.A2007201.h10v05", "MCD43A4.A2002233.h10v05",
           "MCD43A4.A2002241.h10v05")

miss3 <- c("MCD43A4.A2006201.h10v04", "MCD43A4.A2008049.h10v04", "MCD43A4.A2008081.h10v04",
           "MCD43A4.A2012273.h10v04", "MCD43A4.A2012281.h10v04", "MCD43A4.A2012297.h10v04",
           "MCD43A4.A2012321.h10v04", "MCD43A4.A2012353.h10v04", "MCD43A4.A2013001.h10v04")

rescue <- NULL
for(i in seq_along(miss3)){
   rescue[i] <- 
      unique(
         as.character(
            granules[
               which(
                  grepl(
                     pattern = miss3[i],
                     x = granules$Producer.Granule.ID) == 1),
               "Online.Access.URLs"
               ]
         )
      )
}
registerDoParallel(cores=10)
foreach(i = seq_along(rescue)) %dopar%{
  download.file(
    rescue[i],
    destfile = paste(
      "output/",
      substr(x = rescue[i],start = 82,stop = 96),".hdf", sep=""), 
                quiet = TRUE, mode = "wb", method = "wget")
}
# I changed missing file names and directory by hand afetr downloading.

# extracting and mosaicing bands
# b      length      Rad  SNR
# 1 	 620 -  670 	21.8 	128
# 2 	 841 -  876 	24.7 	201
# 3 	 459 -  479 	35.3 	243
# 4 	 545 -  565 	29.0 	228
# 5 	1230 - 1250 	 5.4 	 74
# 6 	1628 - 1652 	 7.3 	275
# 7 	2105 - 2155 	 1.0 	110

# we need bands 2 and 5
# files have same names at each tile
############### #### ### ## # - IMAGE PROCESSING - # ## ### #### ###############
setwd("~/big/USA/MODIS")
library(rgdal)
library(gdalUtils)
library(raster)
library(doParallel)
library(maptools)
############ Create mosaic 
# subset extension in MODIS coordinate system (could be calculated from a shape file)
extent(readShapePoints("K_3_MODproj.shp"))

mask.ext <- c(-8995835.65,4052136.73, -7932725.14,4469948.27) # xmin,ymin , xmax,ymax @QGIS K_3MODproj.shp


f <-list.files("output/h09v05/", pattern = ".hdf$")
# f04 <- list.files("output/h10v04/", pattern = ".hdf$")
# f[which(f %in% f04 == F)]
#pj <- projection(hdfImage[[1]])
rasterOptions(tmpdir="~/big/USA/MODIS/output/temp/")

registerDoParallel(cores=11)
foreach(i = seq_along(f)) %dopar%{
  hdfImage <- 
    list()
  for(j in seq_along(tiles)){
    hdfImage[[j]] <- 
      readGDAL(
        paste0(
          "HDF4_EOS:EOS_GRID:",
          "output/",
          tiles[j],
          "/",
          f[i],
          ":MOD_Grid_BRDF:Nadir_Reflectance_Band5"
          )
        )  
  }
  m <- mosaic(
    raster(hdfImage[[1]]),
    raster(hdfImage[[2]]),
    raster(hdfImage[[3]]),
    fun = mean
  )
  m <- crop(m, mask.ext)
  writeRaster(x = m, filename = paste0("output/bands/",f[i],".B5.tif"), overwrite=TRUE)
  rm(m)
  rm(hdfImage)
  
}
#removeTmpFiles(h=24)
showTmpFiles()
### Estimation NDWI
bands <- c(".B2.tif", ".B5.tif")

foreach(i = seq_along(f)) %dopar%{
  hdfImage <- list()
  for(j in seq_along(bands)) {
    hdfImage[[j]] <- raster(paste0("output/bands/",f[i],bands[j]))    
  }

  ndwi <- (hdfImage[[1]] - hdfImage[[2]])/
          (hdfImage[[1]] + hdfImage[[2]])
  writeRaster(x = ndwi,
              filename = paste0("output/bands/",f[i],".ndwi.tif"),
              overwrite=TRUE)
}

### Estimation of mean and sd by period

# crate list of images from the same period
list.per <- list()
for(i in seq_along(period)){
  list.per[[i]] <- list.files(path = "output/bands/",
                            pattern= paste0(period[i],
                                           ".hdf.ndwi.tif"
                                           )
                            )
}

# stack images from the same period (raster package)
# I found this script here http://markmail.org/thread/lsjnczi5fppddrpr
registerDoParallel(cores=11)
foreach(
  i = seq_along(
    period
  )
) %dopar% {
  s <- 
    stack(
      paste0(
        "output/bands/",
        list.per[[i]]
      )
    )
 # s <- s/10000 # to reduce no significan digits
  # mean and sd  
  m <- 
    mean(
      s,
      na.rm = TRUE
    )
  std <- 
    calc(
      s,
      fun = sd,
      na.rm=TRUE
    )
  # save and delete objects
  writeRaster(
    x = m,
    filename =  
      paste0(
        "output/NDWI/ndwi.",
        period[i],".mean.tif"
      ),
    overwrite=T
  )
  
  writeRaster(
    x = std,
    filename =  
      paste0(
        "output/NDWI/ndwi.",
        period[i],".sd.tif"
      ),
    overwrite = TRUE
  )
  
  rm(s)
  rm(m)
  rm(std)
}

# mean per period
library(raster)
rm(list=ls())
setwd("~/big/USA/MODIS/")
f <-list.files("output/NDWI/", pattern = "mean.tif$")
f[13:15]
# "ndwi.097.mean.tif" "ndwi.105.mean.tif" "ndwi.113.mean.tif"
writeRaster(
  x = mean(
    stack(
      paste0(
        "output/NDWI/",
        f[13:15]
      )
    ),
    na.rm = TRUE
  ),
  filename = 
    "output/NDWI/ndwi.b.spring.tif",
  overwrite = TRUE
)

#f[24:29]
#[1] "ndwi.185.mean.tif" "ndwi.193.mean.tif" "ndwi.201.mean.tif" "ndwi.209.mean.tif" "ndwi.217.mean.tif"
#[6] "ndwi.225.mean.tif"

writeRaster(
  x = mean(
    stack(
      paste0(
        "output/NDWI/",
        f[24:29]
      )
    ),
    na.rm = TRUE
  ),
  filename = 
    "output/NDWI/ndwi.a.summaer.tif",
  overwrite = TRUE
)