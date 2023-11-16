library(RTMB)
dat <- read.table("bh.dat", header = TRUE)
par <- list(loga = 0, logb = 0, logsigma = 0)
nll <- function(par) {
  getAll(par, dat)
  logR <- OBS(logR)
  sigma <- exp(logsigma)
  pred <- loga + log(ssb) - log(1 + exp(logb) * ssb)
  REPORT(pred)
  -sum(dnorm(logR, pred, sigma, TRUE))
}
obj <- MakeADFun(nll, par, silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
odat <- dat
doone <- function() {
  dat$ logR <<- obj$simulate()$ logR
  objs <- MakeADFun(nll, par, silent = TRUE)
  opts <- nlminb(objs$par, objs$fn, objs$gr)
  objs$ report()$ pred
}
sim <- replicate(100, doone())
