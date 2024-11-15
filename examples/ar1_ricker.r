library(RTMB) #install rtmb
data =  read.csv("data/sr_dat.csv")
data$year =  1:nrow(data)
names(data)[2:3] =  c("R", "S") # test

# get starting values 
ls_ricker =  lm(log(data$R) ~ data$S, offset = log(data$S))
ls_report =  summary(ls_ricker)

data =  list(logR = log(data$R), S = data$S)

par =  list(
  log_alpha = coef(ls_ricker)[1],
  eps_a = rep(0, length(data$logR)),
  beta = -coef(ls_ricker)[2],
  log_sd_eps = log(sqrt(ls_report$sigma^2 / 2)),
  log_sd_obs = log(sqrt(ls_report$sigma^2 / 2)), # 1/2 resid var
  trans_rho = 1 # transformed rho
)

# write a function to transform par from - inf to inf -> -1 to 1
to_cor =  function(x) {
  2 / (1 + exp(-2 * x)) - 1
}

f =  function(par) {
  getAll(data, par)
  logR <- OBS(logR)
  n_year <- length(logR)
  rho =  to_cor(trans_rho) # back transform
  sd_eps =  exp(log_sd_eps) # unconditional sd of AR1 process
  sd_obs =  exp(log_sd_obs)
  jnll =  0
  jnll =  jnll - dnorm(eps_a[1], 0, sqrt(1 - rho^2) * sd_eps, TRUE) # initialize
  for (t in 2:n_year) {
    jnll =  jnll - dnorm(eps_a[t],                 # current eps
                         rho * eps_a[t - 1],       # is a function of eps[t-1]
                         sd_eps, # + some stationary noise
                         TRUE
                        )
  }
  pred_log_R =  log_alpha + eps_a - beta * S + log(S)
  alphas = exp(log_alpha + eps_a) 
  ADREPORT(alphas)
  REPORT(rho)
  jnll =  jnll - sum(dnorm(logR, pred_log_R, sd_obs, TRUE))
  jnll 
}

obj =  MakeADFun(f, par, random = c("eps_a"))

opt =  nlminb(obj$par, obj$fn, obj$gr)
opt
sdr =  sdreport(obj)
sdr

obj$simulate()
checkConsistency(obj)

# osa <- oneStepPredict(obj, "logR", "keep")
# qqnorm(osa$residual); abline(0,1)

#----------
# plotting 

year =  1:length(data$logR)
plot(sdr$par.random ~ year,
  ylab = expression(alpha[t]), xlab = "year",
  ylim = c(-3, 3)
)

for (t in 1:length(year)) {
  lines(
    x = rep(t, 2), y = sdr$par.random[t] + c(-1.96, 1.96) * sdr$sd[t],
    col = "#115d9a"
  )
}
