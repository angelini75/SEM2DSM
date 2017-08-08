# test

library(lavaan)
set.seed(1234)
x <- rnorm(n = 50, mean = 20, sd = 3)
y <- rnorm(n = 50, mean = 15, sd = 7)
w <- rnorm(n = 50, mean = 0, sd = 0.5) * x
z <- 0.3 * x + 0.5 * y + 0.5 * w + rnorm(n = 50, mean = 0.2, sd = 2)

data <- data.frame(x = x, y = y, w = w, z = z)

sdn <- function(x){
  n <- length(x)
  m <- sum(x)/n
  var <- sum((x - m)^2)/n
  sd <- sqrt(var)
  sd
}
meann <- function(x){
  n <- length(x)
  m <- sum(x)/n
  m
}

sd.data <- apply(MARGIN = 2, X = data,FUN =  sdn)
mean.data <- apply(MARGIN = 2, X = data,FUN =  meann)

st <- function(x, sd, mean){
  a <- (x - mean)/sd
  a
}

data.s <- data
for(i in seq_along(sd.data)){
  data.s[,i] <- (data[,i] - mean.data[i]) / sd.data[i]
}

data.t <- data
for(i in seq_along(sd.data)){
  data.t[,i] <- (data[,i] - mean(data[,i])) / sd(data[,i])
}



my.model <- '
# X =~ 1 * x
# Y =~ 1 * y
Z =~ 1 * z
z ~~ 0.05 * z

Z ~ w + x + y
'

fit.s <- sem(model = my.model, data = data.s, fixed.x = T)
fit.t <- sem(model = my.model, data = data.t, fixed.x = T)

inspect(object = fit.s, what = "est")
inspect(object = fit.t, what = "est")

lavaan:::computeSigmaHat.LISREL(lavTech(fit.t, "start"))

