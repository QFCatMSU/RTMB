library(RTMB)
data =  read.csv("data/sr_dat.csv")
data$year =  1:nrow(data)
names(data)[1:2] =  c("R", "S")

data =  list(
  S = data$S,
  log_RS_obs = log(data$R / data$S),
  year = data$year
)

pars =  list(
  log_alpha0 = 0,
  beta = 0,
  log_sdr = 0,
  log_sd_alpha = 0,
  alphas = rep(0, length(data$year))
)

f = function(pars) {
  getAll(data, pars)
  jnll =  0 # initialize
  # initial state
  jnll =  jnll - dnorm(alphas[1], exp(log_alpha0), exp(log_sd_alpha), TRUE)
  # walk through the remaining years
  for (t in 2:length(data$year)) {
    jnll =  jnll - dnorm(alphas[t], alphas[t - 1], exp(log_sd_alpha), TRUE)
  }
  # calculate systematic component:
  log_RS_pred =  alphas - beta * S
  # likelihood for observations vs. predictions:
  jnll =  jnll - sum(dnorm(log_RS_obs, log_RS_pred, exp(log_sdr), TRUE))
  REPORT(log_RS_pred)
  ADREPORT(alphas)
  jnll
}

obj =  MakeADFun(f, pars, random = "alphas")
opt =  nlminb(obj$par, obj$fn)
sdr =  sdreport(obj)
opt
sdr 
