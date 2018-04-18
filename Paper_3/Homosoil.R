# HOMOSOIL
# from Chapter 10, R DSM book (Malone, Minasny, McBratney)

rm(list=ls()[])
# install.packages("devtools") 
library(devtools)
# install_bitbucket("brendo1001/ithir/pkg") #ithir package

library(ithir)
library(raster)
library(rasterVis)
data(homosoil_globeDat)


homosoil <- function (grid_data,recipient.lon,recipient.lat) {
  #Gower’s similarity function
  gower <- function(c1, c2, r) 1-(abs(c1-c2)/r)
  #index global data
  grid.lon <-grid_data[,1] #longitude
  grid.lat <-grid_data[,2] #latitude
  grid.climate <-grid_data[,3:54] #climate data
  grid.lith <-grid_data[,58] #lithology
  grid.topo <-grid_data[,55:57] #topography
  #data frame to put outputs of homosoil function
  world.grid<- data.frame(grid_data[, c("X", "Y")],
                          fid = seq(1, nrow(grid_data), 1), homologue = 0, homoclim=0,
                          homolith=0, homotop=0)
  # find the closest recipient point
  dist = sqrt((recipient.lat - grid.lat)^2 + (recipient.lon
                                              - grid.lon)^2)
  imin = which.min(dist)
  
  # climate, lithology and topography for recipient site
  recipient.climate <- grid.climate[imin, ]
  recipient.lith <- grid.lith[imin]
  recipient.topo <- grid.topo[imin, ]
  
  # range of climate variables
  rv <- apply(grid.climate, 2, range)
  rr <- rv[2, ] - rv[1, ]
  # calculate similarity to all variables in the grid
  S <- (mapply(gower, c1 = grid.climate, c2 = recipient.climate,
               r = rr))
  Sr <- apply(S, 1, mean)# take the average
  
  # row index for homoclime with top X% similarity.
  iclim = which(Sr >= quantile(Sr, 0.85), arr.ind = TRUE)
  # save homoclime result
  world.grid$homologue[iclim] <- 1
  world.grid$homoclim[iclim] <- 1
  
  # find within homoclime, areas with homolith
  ilith = which(grid.lith == recipient.lith, arr.ind = TRUE)
  #global comparison
  # homolith in areas of homoclime
  clim.match <- which(world.grid$homologue == 1)
  climlith.match <- clim.match[clim.match %in% ilith]
  
  # save homolith result
  world.grid$homologue[climlith.match] <- 2
  world.grid$homolith[climlith.match] <- 1
  
  # range of topographic variables
  rv <- apply(grid.topo, 2, range)
  rt <- rv[2, ] - rv[1, ]
  # calculate similarity of topographic variables
  Sa <- (mapply(gower, c1 = grid.topo, c2 = recipient.topo, r = rt))
  St <- apply(Sa, 1, mean) # take the average
  # row index for homotop
  itopo = which(St >= quantile(St, 0.85), arr.ind = TRUE)
  
  top.match <- which(world.grid$homologue == 2)
  lithtop.match <- top.match[top.match %in% itopo]
  # save homotop result
  world.grid$homologue[lithtop.match] <- 3
  world.grid$homotop[lithtop.match] <- 1
  
  # homologue raster object
  r1 <- rasterFromXYZ(world.grid[, c(1, 2, 4)])
  r1 <- as.factor(r1)
  rat <- levels(r1)[[1]]
  rat[["homologue"]] <- c("", "homoclime", "homolith", "homotop")
  levels(r1) <- rat
  
  retval <- list(r1, world.grid)
  return(retval)}


recipient.lat = -33.971334
recipient.lon = -60.021632
xy1 <- structure(list(longitude = recipient.lon, latitude = recipient.lat), 
                  .Names = c("longitude", "latitude"), 
                  class = "data.frame", row.names = c(NA, -1L))

xy2 <- structure(list(longitude = -98.621095, latitude = 39.175715), 
                .Names = c("longitude", "latitude"), 
                class = "data.frame", row.names = c(NA, -1L))

ba <- SpatialPointsDataFrame(coords = xy1, data = xy1,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

ks <- SpatialPointsDataFrame(coords = xy2, data = xy2,
                             proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))


result <- homosoil(grid_data = homosoil_globeDat,
                   recipient.lon = recipient.lon,
                   recipient.lat = recipient.lat)


# plot
area_colors <- c("#EFEFEF", "#666666", "#FFDAD4", "#FF0000")
levelplot(result[[1]], col.regions = area_colors,
          xlab = "", ylab = "") + 
  layer(sp.points(ba, col = "green", pch = 20, cex = 2))+ 
  layer(sp.points(ks, col = "green", pch = 20, cex = 2))



  