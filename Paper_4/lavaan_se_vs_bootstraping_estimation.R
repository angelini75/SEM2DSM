# test bootstraping to test how robust are the coefficients
rm(list=ls()[])
library(lavaan)
library(doParallel)
library(reshape)
library(ggplot2)

load("~/Documents/SEM2DSM1/Paper_4/data/se_vs_bootstraping.RData")
ks <- read.csv("~/Documents/SEM2DSM1/Paper_3/data/KS.data-0.3.csv")[,-1: -2]
names(ks)[1:9] <- c("Clay.A","Clay.B","Clay.C","CEC.A","CEC.B","CEC.C","OC.A","OC.B","OC.C")
s <- ks
# get the estimates  
registerDoParallel(cores = 4)
trials = 1000
x <- foreach(icount(trials), .combine=rbind) %dopar% {
  samples <- s[sample(x = 1:nrow(s), size = nrow(s), replace = TRUE),]
  rownames(samples) <- 1:147
  fit <- sem(model = my.model.lv,data = samples, meanstructure = FALSE, 
             fixed.x = T)
  lavaan:::lav_model_get_parameters(lavmodel = fit@Model)
}
x <- as.data.frame(x)
# estimates of the model (fit) published in paper 3
my.fit <- sem(model = my.model.lv,data = ks, meanstructure = FALSE, 
              fixed.x = T)

# names of the estimates
p <- partable(my.fit)
colnames(x) <- paste0(p$lhs[p$free!=0],p$op[p$free!=0],p$rhs[p$free!=0])
# 

# statistics of the estimates
xmean <- colMeans(x)
xsd <- apply(x, 2, sd)

xmean <- reshape::melt(as.data.frame(xmean))
xmean$variable <- names(x)

my.x <- lavaan:::lav_model_get_parameters(lavmodel = my.fit@Model)
my.x <- reshape::melt(as.data.frame(my.x))
my.x$variable <- names(x)

# compare bootstrap mean estimates with my.fit estimates
ggplot2::ggplot(data = reshape::melt(x), mapping = aes(x = value)) +
  geom_histogram(bins = 40) + facet_wrap(~variable, scales = 'free_x') +
  geom_vline(aes(xintercept=value, color='my.x'), data = my.x) +
  geom_vline(aes(xintercept=value, color='bootstrap'), data=xmean)

p <- partable(my.fit)
p$se.boots <- NA
p$se.boots[p$free != 0] <- xsd
p$se[p$free != 0] - p$se.boots[p$free != 0]

lav_model_x2GLIST <- function (lavmodel = NULL, x = NULL) 
{
  GLIST <- lavmodel@GLIST
  for (mm in 1:length(GLIST)) {
    M.EL.IDX <- lavmodel@m.free.idx[[mm]]
    X.EL.IDX <- lavmodel@x.free.idx[[mm]]
    GLIST[[mm]][M.EL.IDX] <- x[X.EL.IDX]
  }
  GLIST
}
fit2 <- lav_model_x2GLIST(lavmodel = fit@Model,x = x)

fit@Model@GLIST$beta[1:9,1:9]
fit2$beta[1:9,1:9]

