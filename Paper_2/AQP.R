rm(list=ls())
name <- function(x) { as.data.frame(names(x))}
# chose one
#setwd("~/big/SEM2DSM1/Paper_2/data/")
setwd("~/Documents/SEM2DSM1/Paper_2/data/")

d <- read.csv("calib.data-sp.csv")
p <- unique(d)
p$top[p$id.p.h=="350_B1"]  <- 25
p[p$id.p==568,1:10][2,4] <- 18
p[p$id.p==576,1:10][5,4] <- 42
p[p$id.p==710,1:10][7,c(4,5)] <- c(160,180)
p[p$id.p==687,1:10][c(5,6),c(4,5)] <- c(145,100,190,145)
p$hor[p$hor== "AB|BA"] <- "AB.BA"
p$name <- as.factor(p$hor)
#p <- d[d$id.p== c(619, 350),1:15]

library(aqp)
library(lattice)
aqp::depths(p)<- id.p ~ top + bottom 

name(p)
p[3,]
horizons(p)

a <- slab(p, fm = ~ CEC)
b <- slab(p, fm = ~ OC )
c <- slab(p, fm = ~ clay )
# stack into long format
abc <- lattice::make.groups(a, b, c)
abc$which <- factor(abc$which, levels=c('a','b', "c"), 
                   labels=c('CEC', 'OC', "Clay"))

setwd("~/Dropbox/PhD Marcos/Paper 2/Figures/")

tiff(filename = "outfile.tif", width = 1600, height = 1200, res =  300)
lattice::xyplot(top ~ p.q50 | which, data=abc, ylab='Depth (cm)',
       xlab='Median bounded by 25th and 75th percentiles',
       lower=abc$p.q25, upper=abc$p.q75, ylim=c(250,-5),
       panel=panel.depth_function, 
       prepanel=prepanel.depth_function,
       cf=abc$contributing_fraction,
       layout=c(3,1), scales=list(x=list(alternating=1, relation="free")))
dev.off() 
# select a profile to use as the basis for simulation
s <- p[3, ]

# reset horizon names
s$name <- paste('H', seq_along(s$hor), sep='')

# simulate 25 new profiles, using 's' and function defaults
sim.1 <- sim(s, n=25)

# simulate 25 new profiles using 's' and variable SD for each horizon
sim.2 <- sim(s, n=25, hz.sd=c(1, 2, 5, 5, 5, 2))

par(mfrow=c(2,1), mar=c(0, 0, 0, 0))
plot(sim.1)
mtext('SD = 2', side=2, line=-1.5, font=2, cex=0.75)
plot(sim.2)
mtext('SD = c(1, 2, 5, 5, 5, 2)', side=2, line=-1.5, font=2, cex=0.75)

# aggregate horizonation of simulated data
# note: set class_prob_mode=2 as profiles were not defined to a constant depth
# sim.2$name <- factor(sim.2$name)
# p$name <- factor(p$hor)
a <- slab(p[p$env==3], ~ name, class_prob_mode=2)

# convert to long format for plotting simplicity
library(reshape)
a.long <- melt(a, id.vars=c('top','bottom'), measure.vars=levels(p$name))

# plot horizon probabilities derived from simulated data
# dashed lines are the original horizon boundaries
library(lattice)
xyplot(top ~ value, groups=variable, data=a.long, subset=value > 0,
       ylim=c(200, -5), type=c('l','g'), asp=1.5,
       ylab='Depth (cm)', xlab='Probability', 
       auto.key=list(columns=4, lines=TRUE, points=FALSE),
       panel=function(...) {
         panel.xyplot(...)
         panel.abline(h=s$top, lty=2, lwd=2)
       })


sim.1@horizons


