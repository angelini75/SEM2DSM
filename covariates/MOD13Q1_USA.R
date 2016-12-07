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

# functions to read directories and files
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
MODISP <- "MOD13Q1"
# driver and bands of HDF file for GDAL
driver <- "HDF4_EOS:EOS_GRID:"
bands <- data.frame(name = c(":MODIS_Grid_16DAY_250m_500m_VI:250m 16 days EVI"),
                    file.name = c("EVI"))
#define tiles
tiles <- c("h09v05", "h10v04", "h10v05")
#define years
yrs <- as.character(2000:2015)
#get periods
library(RCurl)
period <- gsub("ftp://ladsweb.nascom.nasa.gov/allData/5/MOD13Q1/2005/","",
               readdir("ftp://ladsweb.nascom.nasa.gov/allData/5/MOD13Q1/2005/"))

readdir("152.61.133.130/MOLT/MOD13Q1.005")

# Files to be download ####
# http://reverb.echo.nasa.gov/reverb/
MODProd <- read.csv(file = "MOD13Q1_reverb.csv")

# two dataframes, one for each tile 
h09v05 <- MODProd[grep(pattern = "?h09v05",x = MODProd$Producer.Granule.ID),]
h10v04 <- MODProd[grep(pattern = "?h10v04",x = MODProd$Producer.Granule.ID),]
h10v05 <- MODProd[grep(pattern = "?h10v05",x = MODProd$Producer.Granule.ID),]
h09v05$Online.Access.URLs <- as.character(h09v05$Online.Access.URLs)
h10v04$Online.Access.URLs <- as.character(h10v04$Online.Access.URLs)
h10v05$Online.Access.URLs <- as.character(h10v05$Online.Access.URLs)

# Download the files ####
registerDoParallel(cores=6)
foreach(i = seq_along(h09v05[,1])) %dopar%{
  download.file(h09v05[i,5], 
                destfile = paste("output/",
                                 tiles[1],
                                 "/",
                                 substr(x = h09v05[i,2],start = 9,stop = 16),
                                 ".hdf",
                                 sep=""),
                quiet = TRUE, 
                mode = "wb",
                method = "wget",
                extra = "--load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies")
}
foreach(i = seq_along(h10v04[,1])) %dopar%{
  download.file(h10v04[i,5], 
                destfile = paste("output/",
                                 tiles[2],
                                 "/",
                                 substr(x = h10v04[i,2],start = 9,stop = 16),
                                 ".hdf",
                                 sep=""),
                quiet = TRUE, 
                mode = "wb",
                method = "wget",
                extra = "--load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies")
}
foreach(i = seq_along(h10v05[,1])) %dopar%{
  download.file(h10v05[i,5], 
                destfile = paste("output/",
                                 tiles[3],
                                 "/",
                                 substr(x = h10v05[i,2],start = 9,stop = 16),
                                 ".hdf",
                                 sep=""),
                quiet = TRUE, 
                mode = "wb",
                method = "wget",
                extra = "--load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies")
}

############### #### ### ## # - IMAGE PROCESSING - # ## ### #### ###############
setwd("~/big/USA/MODIS")
library(rgdal)
library(gdalUtils)
library(raster)
library(doParallel)
library(maptools)
############ Create mosaic 
# subset extension in MODIS coordinate system (could be calculated from a shape file)
aoi <- readShapePoly("Platte_area.shp")

# define projections
wgs84 <- CRS("+init=epsg:4326")
modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# assign projection
proj4string(aoi) <- wgs84

# use spTransform
aoi <- spTransform(aoi, modis)



f <-list.files("output/h10v04/", pattern = ".hdf$")
# f04 <- list.files("output/h10v04/", pattern = ".hdf$")
# f[which(f %in% f04 == F)]
#pj <- projection(hdfImage[[1]])
rasterOptions(tmpdir="~/big/USA/MODIS/output/temp/")

registerDoParallel(cores=11)
foreach(i = seq_along(f)) %dopar%{
  for(k in seq_along(bands[,1])){
    hdfImage <- list()
    for(j in seq_along(tiles)){
      hdfImage[[j]] <- readGDAL(
        paste0("HDF4_EOS:EOS_GRID:",
               "output/",
               tiles[j],
               "/",
               f[i],
               bands[k,1]
        )
      )  
    }
    m <- mosaic(
      raster(hdfImage[[1]]),
      raster(hdfImage[[2]]),
      fun = mean
    )
    m <- crop(m, aoi)
    writeRaster(x = m, 
                filename = paste0("output/bands/",
                                  f[i],
                                  ".",
                                  bands[k,2],
                                  ".tif"),
                overwrite=TRUE)
    rm(m)
    rm(hdfImage)
  }
}
#removeTmpFiles(h=24)
showTmpFiles()

### Estimation of mean and sd by period
# crate list of images from the same period
list.per <- list()
for(i in seq_along(period)){
  list.per[[i]] <- list.files(path = "output/bands/",
                              pattern= paste0(period[i],
                                              ".hdf.EVI.tif"
                              )
  )
}

# stack images from the same period (raster package)
# I found this script here http://markmail.org/thread/lsjnczi5fppddrpr
registerDoParallel(cores=11)
foreach(i = seq_along(period)) %dopar% {
  s <- stack(paste0("output/bands/",
                    list.per[[i]]))
  # s <- s/10000 # to reduce no significan digits
  # mean and sd  
  m <- mean(s, na.rm = TRUE)
  std <- calc(s, fun = sd, na.rm=TRUE)
  
  # save and delete objects
  writeRaster(x = m,
              filename = paste0("output/EVI/EVI.",
                                period[i],
                                ".mean.tif"),
              overwrite=T)
  
  writeRaster(x = std,
              filename =  paste0("output/EVI/EVI.",
                                 period[i],
                                 ".sd.tif"),
              overwrite = TRUE)
  
  rm(s)
  rm(m)
  rm(std)
}

# crate list of images from entire period
EVI <- list.files(path = "output/bands/",
                  pattern= ".hdf.EVI.tif")

s <- stack(paste0("output/bands/", EVI))
# mean and sd  
m <- mean(s, na.rm = TRUE)
std <- calc(s, fun = sd, na.rm=TRUE)

# save and delete objects
writeRaster(x = m, filename = "output/EVI/EVI.mean.tif", overwrite=TRUE)

writeRaster(x = std, filename = "output/EVI/EVI.sd.tif", overwrite = TRUE)

