rm(list=ls())
name <- function(x) { as.data.frame(names(x))}
panel.depth_function_MA <- 
  function (x, y, id, upper = NA, lower = NA, subscripts = NULL, 
            groups = NULL, sync.colors = FALSE, cf = NA, cf.col = NA, 
            cf.interval = 20, ...) 
  {
    panel.grid(h = -1, v = -1, lty = 3, col = "#999999")
    superpose.line <- trellis.par.get("superpose.line")
    if (length(y) > length(x)) {
      if (missing(id)) 
        stop("must provide a profile id")
      message("plotting segments...")
      if (!missing(groups)) 
        d <- data.frame(prop = x, bnd = y, upper = upper[subscripts], 
                        lower = lower[subscripts], groups = groups[subscripts], 
                        id = id[subscripts])
      else d <- data.frame(prop = x, bnd = y, upper = upper[subscripts], 
                           lower = lower[subscripts], groups = factor(1), id = id[subscripts])
      by(d, d$id, make.segments, ...)
    }
    else {
      if (!missing(upper) & !missing(lower)) {
        if (!missing(groups) & !missing(subscripts)) {
          d <- data.frame(yhat = x, top = y, upper = upper[subscripts], 
                          lower = lower[subscripts], groups = groups[subscripts])
          ll <- levels(d$groups)
          n_groups <- length(ll)
        }
        if (missing(groups)) {
          d <- data.frame(yhat = x, top = y, upper = upper[subscripts], 
                          lower = lower[subscripts], groups = factor(1))
          ll <- levels(d$groups)
          n_groups <- length(ll)
        }
        if (sync.colors) 
          region.col <- rep(superpose.line$col, length.out = n_groups)
        else region.col <- rep(grey(0.7), length.out = n_groups)
        by(d, d$groups, function(d_i) {
          m <- match(unique(d_i$group), ll)
          d_i <- d_i[which(!is.na(d_i$upper) & !is.na(d_i$lower)), 
                     ]
          panel.polygon(x = c(d_i$lower, rev(d_i$upper)), 
                        y = c(d_i$top, rev(d_i$top)), col = region.col[m], 
                        border = NA, ...)
        })
      }
      else {
        d <- data.frame(yhat = x, top = y, groups = groups[subscripts])
        ll <- levels(d$groups)
        n_groups <- length(ll)
      }
      line.col <- rep(superpose.line$col, length.out = n_groups)
      line.lty <- rep(superpose.line$lty, length.out = n_groups)
      line.lwd <- rep(superpose.line$lwd, length.out = n_groups)
      by(d, d$groups, function(d_i) {
        m <- match(unique(d_i$group), ll)
        panel.lines(d_i$yhat, d_i$top, lwd = line.lwd[m], 
                    col = line.col[m], lty = line.lty[m])
      })
    }
    if (!missing(cf)) {
      d$cf <- cf[subscripts]
      by(d, d$groups, function(d_i) {
        m <- match(unique(d_i$group), ll)
        if (is.na(cf.col)) {
          cf.col <- line.col[m]
        }
        cf.approx.fun <- approxfun(d_i$top, d_i$cf, method = "linear")
        y.q95 <- quantile(d_i$top, probs = c(0.95), na.rm = TRUE)
        a.seq <- seq(from = 2, to = y.q95, by = cf.interval)
        a.seq <- a.seq + ((m - 1) * cf.interval/4)
        a.CF <- cf.approx.fun(a.seq)
        a.text <- paste(round(a.CF * 100), "%")
        not.na.idx <- which(!is.na(a.CF))
        a.seq <- a.seq[not.na.idx]
        a.text <- a.text[not.na.idx]
        unit <- gpar <- NULL
      })
    }
  }

# chose one
#setwd("~/big/SEM2DSM1/Paper_2/data/")
setwd("~/Documents/SEM2DSM1/Paper_2/data/")

d <- read.csv("calib.data-sp.csv")
p <- unique(d)
p[p$id.p==568,1:10][2,4] <- 18
p[p$id.p==576,1:10][5,4] <- 42
p <- p[p$hor!= "BC",]
p$name <- as.factor(as.character(p$hor))

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


sp <- xyplot(top ~ p.q50 | which, data=abc, ylab = "", xlab = "",  #ylab='Depth (cm)',
       #xlab='Median bounded by 25th and 75th percentiles',
       lower=abc$p.q25, upper=abc$p.q75, ylim=c(250,-5), type=c('l', "g"),
       panel=panel.depth_function_MA,
       prepanel=prepanel.depth_function,
       cf=abc$contributing_fraction,
       layout=c(3,1), scales=list(x=list(alternating=1, relation="free")),
       strip = strip.custom(
         factor.levels = c(expression("CEC"~~"/ cmol"[c]~~"kg"^{-1}),
                           "OC / %",
                           "Clay / %")),
       par.settings=list(grid.pars=list(fontfamily="serif")),
       sub = expression("Median bounded by 25"^{th}~"and 75"^{th}~"percentiles"))


a <- slab(p, ~ name, class_prob_mode=2)
s <- p[3, ]
# convert to long format for plotting simplicity
library(reshape)
a.long <- melt(a, id.vars=c('top','bottom'), measure.vars=levels(p$name))

# plot horizon probabilities derived from simulated data
# dashed lines are the original horizon boundaries
library(lattice)
hor <- xyplot(top ~ value, groups=variable, data=a.long, subset=value > 0,
       ylim=c(250, -5), type=c('l', "g"), asp=1.5,
       ylab='', xlab='',draw.key = TRUE,
       #panel=panel.depth_function,
       panel=panel.depth_function_MA,
       cf=abc$contributing_fraction, cex=0.5, 
       strip = strip.custom(which.given = 3,
                            strip.names = TRUE,
                            var.name = "Horizons"),
       par.settings=list(grid.pars=list(fontfamily="serif")),
       auto.key=list(columns=3, lines=TRUE,
                     points = FALSE, space = "bottom"))
tiff(filename = "~/Dropbox/PhD Marcos/Paper 2/Figures/sp.tif",
     width = 1500, height = 1300, res =  300)
sp
dev.off()
tiff(filename = "~/Dropbox/PhD Marcos/Paper 2/Figures/hor.tif",
     width = 1500, height = 1300, res =  300)
hor
dev.off()

library(latticeExtra)

tiff(filename = "~/Dropbox/PhD Marcos/Paper 2/Figures/Fig2.tif",
     width = 2500, height = 1300, res =  320)
update(c(sp,hor),
       xlab = NULL, ylab = "Depth / cm", x.same = FALSE, 
       y.same = TRUE, layout = c(4, 1), 
         scales = list(alternating= FALSE),
         strip = strip.custom(
           factor.levels = c(expression("CEC"~~"/ cmol"[c]~~"kg"^{-1}),
                             "OC / %",
                             "Clay / %",
                             "Horizons / frequency")))
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


sim.1@horizons


