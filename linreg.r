# install.packages("RTMB", repos = c("https://kaskr.r-universe.dev", "https://cloud.r-project.org"))
library(RTMB)

#------------------------------------------
# simulate and estimate a linear regression
n = 100 # n data
sd = 5 # sd of likelihood
b0 = 10 # intercept
b1 = 5 # slope

set.seed(123)
x1 = runif(n, -1, 1)
y_obs = rnorm(n, b0 + b1 * x1, sd)
#-------------------------------------------
# estimate it

data = list(
  n = n,
  y_obs = y_obs,
  x1 = x1
)

pars = list(
  b0 = 1,
  b1 = 1,
  log_sd = log(3)
)

f = function(pars) {
  getAll(data, pars) # replaces DATA_XX, PARAMETER_YY
  y_pred = b0 + b1 * x1
  nll = -sum(dnorm(y_obs, y_pred, exp(log_sd), log = TRUE))
  nll
}

obj = MakeADFun(f, pars)
obj$fn() # objective function
obj$gr() # gradients

opt = nlminb(obj$par, obj$fn, obj$gr)
opt

sdr = sdreport(obj)
sdr
