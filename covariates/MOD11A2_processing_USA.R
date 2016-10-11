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
MODISP <- "MOD11A2"
# driver and bands of HDF file for GDAL
driver <- "HDF4_EOS:EOS_GRID:"
band1 <- ":MODIS_Grid_8Day_1km_LST:LST_Day_1km"
#band2 <- ":MOD_Grid_BRDF:Nadir_Reflectance_Band5"
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
    u <-append(u,readdir(paste0(URL,"/",MODISP,"/", yrs[i],"/")))
}
u <- paste0(u,"/")

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

u1
u2
u3

# rescuing missing files
granules <- read.csv(file = "MOD11A2_2001-02-05-06-07-08-09-11-13-14.csv")
miss <- c("MOD11A2.A2001353.h09v05", "MOD11A2.A2002001.h09v05", "MOD11A2.A2002017.h09v05",
          "MOD11A2.A2002145.h09v05", "MOD11A2.A2002153.h09v05", "MOD11A2.A2005353.h09v05",
          "MOD11A2.A2006033.h09v05", "MOD11A2.A2009289.h09v05", "MOD11A2.A2009321.h09v05",
          "MOD11A2.A2001225.h10v04", "MOD11A2.A2001233.h10v04", "MOD11A2.A2001249.h10v04",
          "MOD11A2.A2001281.h10v04", "MOD11A2.A2001289.h10v04", "MOD11A2.A2001345.h10v04",
          "MOD11A2.A2005017.h10v04", "MOD11A2.A2005025.h10v04", "MOD11A2.A2008281.h10v04",
          "MOD11A2.A2009337.h10v04", "MOD11A2.A2011193.h10v04", "MOD11A2.A2011225.h10v04",
          "MOD11A2.A2011233.h10v04", "MOD11A2.A2014249.h10v04", 
          "MOD11A2.A2007281.h10v05", "MOD11A2.A2007289.h10v05", "MOD11A2.A2013233.h10v05",
          "MOD11A2.A2013353.h10v05")

# "MOD11A2.A2001177.h10v05" cannot be rescued so A2001177 is deleted from all tiles

granules$Online.Access.URLs <- as.character(granules$Online.Access.URLs)
rescue <- NULL
for(i in seq_along(miss)){
   rescue[i] <- 
      unique(
         as.character(
            granules[
               which(
                  grepl(
                     pattern = miss[i],
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
#pj <- projection(hdfImage[[1]])

registerDoParallel(cores=10)
foreach(i = seq_along(f)) %dopar%{
  hdfImage <- list()
  for(j in seq_along(tiles)){
    hdfImage[[j]] <- readGDAL(paste0(driver,
                                     "output/",
                                     tiles[j],
                                     "/",
                                     f[i],
                                     band1))  
  }
  m <- mosaic(
     raster(hdfImage[[1]]),
     raster(hdfImage[[2]]),
     raster(hdfImage[[3]]),
     fun = mean
     )
  m <- crop(m, mask.ext)
  writeRaster(x = m, filename = paste0("output/bands/",f[i],".LST.tif"), overwrite=TRUE)
}

### Estimation of mean and sd by period

# crate list of images from the same period
list.per <- list()
for(i in seq_along(period)){
  list.per[[i]] <- 
    list.files(
      path = "output/bands/",
      pattern= paste0(
        period[i],
        ".hdf.LST.tif"
        )
    )
}

# stack images from the same period (raster package)
# I found this script here http://markmail.org/thread/lsjnczi5fppddrpr
registerDoParallel(cores=12)
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
        "output/LST/lst.",
        period[i],".mean.tif"
      ),
    overwrite=T
  )
  
  writeRaster(
    x = std,
    filename =  
      paste0(
        "output/LST/lst.",
        period[i],".sd.tif"
      ),
    overwrite = TRUE
  )
  
  rm(s)
  rm(m)
  rm(std)
}

# get lstm (same as paper 2)
library(raster)
rm(list=ls())
setwd("~/big/USA/MODIS/")
f <-list.files("output/LST/", pattern = "mean.tif$")
f

writeRaster(
  x = mean(
    stack(
      paste0(
        "output/LST/",
        f
      )
    ),
    na.rm = TRUE
  ),
  filename = 
    "output/LST/lstm.tif",
  overwrite = TRUE
)

