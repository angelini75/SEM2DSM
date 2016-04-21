############### #### ### ## # -  - # ## ### #### ###############
# Purpose        : Replace NAs and restructurate calibration data
# Maintainer     : Marcos Angelini  (marcos.angelini@wur.nl); 
# Contributions  : 
# Status         : 
# Note           : result of reshape_data_all_hor.0.2.R
# sessionInfo(@RStudio desktop)  lenovo ThinkPad T430 (4 cores)
# R version 3.0.2 (2013-09-25)
# Platform: x86_64-pc-linux-gnu (64-bit)
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base   
# other attached packages:
#   [1] questionr_0.5
library(questionr)

rm(list=ls())
name <- function(x) { as.data.frame(names(x))}
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration")
sp <- read.csv("calib.data-sp.csv")[,-1]
covar <- read.csv("calib.data-covar.csv")[,-1]

########### Dictionary ###########
# covar = covariates extracted from calib with id, X and Y coordinates
# sp = soil properties (all)
# n 
# x 
#################################
# function to compute hz-thickness weighted mean of Clay, OC and CEC
wt.mean.properties <- function(data, properties) {
  # use horizon thickness as a weight
  d = data
  p = properties
  site <- unique(d$id.p)  
  hor <- c("A", "B", "C")
  m1 <- as.data.frame(matrix(data = NA,nrow = 0,ncol = length(p)+1,
                             dimnames = list(NULL,c("id.p",paste(p[1],rep(hor,1), sep=".")))))
  m2 <- as.data.frame(matrix(data = NA,nrow = 0,ncol = length(p)+1,
                             dimnames = list(NULL,c("id.p",paste(p[2],rep(hor,1), sep=".")))))
  m3 <- as.data.frame(matrix(data = NA,nrow = 0,ncol = length(p)+1,
                             dimnames = list(NULL,c("id.p",paste(p[3],rep(hor,1), sep=".")))))
  for(i in seq_along(site)){
    d.sub <- d[d$id.p == site[i],]
    for(j in seq_along(hor)){
      subset <- d.sub[d.sub$hor == hor[j],]
      thick <- subset$bottom - subset$top  
      # compute the weighted mean, accounting for the possibility of missing data
      m1[i,1] <- site[i]
      m2[i,1] <- site[i]
      m3[i,1] <- site[i]
      m1[i,j+1] <- wtd.mean(subset[which(names(d) == p[1])], weights=thick, na.rm=TRUE)
      m2[i,j+1] <- wtd.mean(subset[which(names(d) == p[2])], weights=thick, na.rm=TRUE)
      m3[i,j+1] <- wtd.mean(subset[which(names(d) == p[3])], weights=thick, na.rm=TRUE)
    }
  }
  r <- merge(merge(m1, m2, by = "id.p", all = TRUE), m3, by = "id.p", all = TRUE)
  r
}
#m <- wt.mean.properties(data = s, properties = c("CEC", "OC","clay"))


