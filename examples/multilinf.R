# Bence July 2023, cahill added a bit of styling...
# R script for fitting vonB model to multiple ponds with TMB
# Allowing for random effects

library(RTMB)

mvonBRdat = read.table("data/multiLinfTMB.dat", header = F)

# Set up the data and starting value of parameters for TMB
data = list(L = as.matrix(mvonBRdat), A = 2:12)

Linfs = sapply(mvonBRdat, max)
pars = list(
  logLinfmn = 7, loglogLinfsd = -2.8,
  logLinfs = log(Linfs),
  logK = -1.6, t0 = 0, logSig = 4
)

f = function(pars) {
  getAll(data, pars)
  Linfmn = exp(logLinfmn)
  logLinfsd = exp(loglogLinfsd)
  Linfs = exp(logLinfs)
  K = exp(logK)
  Sig = exp(logSig)
  nponds = length(Linfs)
  nages = length(A)
  predL = matrix(0, nrow = nages, ncol = nponds)
  # fill one column (pond) at a time:
  for (i in 1:nponds) {
    predL[, i] = Linfs[i] * (1 - exp(-K * (A - t0)))
  }
  nll = -sum(dnorm(x = L, mean = predL, sd = Sig, TRUE))
  nprand = -sum(dnorm(x = logLinfs, mean = logLinfmn, sd = logLinfsd, TRUE))
  jnll = nll + nprand
  jnll
}

obj = MakeADFun(f, pars, random = c("logLinfs"))
fit = nlminb(obj$par, obj$fn, obj$gr)

sdr = sdreport(obj)
sdr