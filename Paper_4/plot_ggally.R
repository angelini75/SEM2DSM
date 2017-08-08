rm(list=ls())
name <- function(x) { as.data.frame(names(x))}
# chose one
# Libraries ####
library(lavaan)
library(pastecs)
library(utils)

setwd("~/Documents/SEM2DSM1/Paper_3/data/")
d <- read.csv("KS.data-0.2.csv")[,c(-1)] 
name(d)
names(d)[5:10] <- c("CEC.A","CEC.B","CEC.C","OC.A","OC.B","OC.C")
# remove outlayers
# d <- d[d$idp!=26058,]
# d <- d[d$idp!=22961,]
d <- cbind(d[1],
           d[,colnames(E)[2:10]],
           d[11:21])
# Descriptive statistics and normality test. ####
round(stat.desc(d,norm = TRUE),3)
# Soil properties does not present strong deviation from normality.
# But some covariates need to be transformed. First, we store original mean and 
# sd in ST
ST <- t(stat.desc(d,norm = TRUE)[c(9,13),])

# Based on normtest.W the following covariates need to be transformed
d$twi <- log10(d$twi)
d$vdchn <- log10(d$vdchn+10)
d$ndwi.a <- (d$ndwi.a+10)^.3
# OC as log10 of OC
# d$OC.A <- log10(d$OC.A)
# d$OC.B <- log10(d$OC.B)
# d$OC.C <- log10(d$OC.C)

# New statistics
round(stat.desc(d,norm = TRUE),3)
# New mean and sd
STt <- t(stat.desc(d,norm = TRUE)[c(9,13),])

# standardised data set ####
std <- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] - st[i,1]) / st[i,2]
  }
  y
}
D <- std(d,STt)
D[,1] <- d[,1] 



library(reshape)
library(ggplot2)
meltp <- unique(melt(d,id.vars = c("idp")))

ggplot(data = meltp,
       aes(x = value, fill=variable)) + geom_histogram() +
  facet_wrap( ~ variable, scales = "free_x")
d$H <- NA
name(d)
e <- data.frame(d[,c(1,2,5,8,11:19)])
e$H <- "A"
name(e)
names(e)[2:4] <- c("clay.B","CEC.B","OC.B")
e <- rbind(e,d[,c(1,3,6,9,11:20)])
name(e)
e$H[is.na(e$H)] <- "B"

names(e)[2:4] <- c("clay.C","CEC.C","OC.C")
e <- rbind(e,d[,c(1,4,7,10,11:20)])
e$H[is.na(e$H)] <- "C"
names(e)[2:4] <- c("Clay","CEC","OC")

meltp <- melt(e,id.vars = c("idp","H"))
ggplot(data = meltp[meltp$variable=="CEC" |
                      meltp$variable=="OC" |
                      meltp$variable=="Clay",],
       aes(x = value, fill = H)) + geom_density(alpha = 0.4) +
  facet_wrap( ~ variable,scales = "free")
ggplot(data = unique(meltp[!(meltp$variable=="CEC" |
                               meltp$variable=="OC" |
                               meltp$variable=="Clay"),]),
       aes(x = value)) + geom_density(alpha = 0.4) +
  facet_wrap( ~ variable,scales = "free")
#############

# # New statistics
# d <- d[,-20:-21]
# names(d)[12:13] <- c("twi", "vdchn")
# round(stat.desc(d,norm = TRUE),3)
# # New mean and sd
# STt <- t(stat.desc(d,norm = TRUE)[c(9,13),])
# 
# # standardised data set ####
# std <- function(x, st){
#   y <- x
#   for(i in seq_along(names(x))){
#     y[,i] <- (x[,i] - st[i,1]) / st[i,2]
#   }
#   y
# }
# 
# D <- std(d,STt)
# 
# D[,1] <- d[,1] 
# # statistics
# round(stat.desc(D,norm = TRUE),0)

#### correlogram
library(psych)
corr <- d[,1:10]
par(mfrow=c(1,1),pty = "s")
pairs.panels(d[2:10], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)#, cex= 1)  # plot correlogram
# ## correlogram ####
library(GGally)
library(ggplot2)
d <- d[d$OC.C<5,]
d <- d[,c(2:13,15:16,18:21)]

png(filename =  "~/Documents/SEM2DSM1/Paper_4/data/plot.1.png",
    width = 5000,height = 5000, res = 400 )
ggscatmat(d, columns = 1:18,  alpha=0.15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()
