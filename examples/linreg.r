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
  y_obs = y_obs,
  x1 = x1
)

par = list(
  b0 = 1,
  b1 = 1,
  log_sd = log(3)
)

f = function(par) {
  getAll(data, par) 
  y_obs <- OBS(y_obs)
  sigma = exp(log_sd)
  y_pred = b0 + b1 * x1
  REPORT(y_pred)
  -sum(dnorm(y_obs, y_pred, sigma, TRUE))
}

obj = MakeADFun(f, par)
obj$fn() # objective function
obj$gr() # gradients

opt = nlminb(obj$par, obj$fn, obj$gr)
opt

sdr = sdreport(obj)
sdr 

# do a simulation 

obj$simulate()

odat <- data 

doone <- function() {
  data$y_obs <<- obj$simulate()$y_obs
  objs <- MakeADFun(f, par, silent = TRUE)
  opts <- nlminb(objs$par, objs$fn, objs$gr)
  objs$report()$y_pred
}

sim <- replicate (100 , doone () )

