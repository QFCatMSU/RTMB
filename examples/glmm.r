library(RTMB)

#-------------------------------------------------------------------
# simulate and estimate a Poisson GLMM
# y_i ~ Poisson(lambda_i)
# log(lambda_i) = b0 + b1x1 + eps_site
# where eps_site ~ N(0, sd_eps)

#--------------------------------------------------------------------

b0 = 3
b1 = 2.3

sd_eps = 0.7
n_site = 15
site = rep(1:n_site, each = n_site)
nobs = length(site)

set.seed(333)
x1 = runif(nobs, -1, 1)
eps_i = rnorm(n_site, 0, sd_eps)
lambda_i = rep(NA, length(site))
nobs = length(lambda_i)
for(i in 1:nobs){
  lambda_i[i] = exp(b0 + b1 * x1[i] + eps_i[site[i]])
}
yobs = rpois(nobs, lambda_i)

#-------------------------------------------------------------------
# estimate it

data = list(
  yobs = yobs,
  x1 = x1, 
  site = site
)

pars = list(
  b0 = 1,
  b1 = 1,
  log_sd_site = log(1),
  eps_site = rep(0, n_site)
)

f = function(pars) {
  getAll(data, pars)
  sd_site = exp(log_sd_site)                           # back transform
  jnll = 0                                             # initialize
  jnll = jnll - sum(dnorm(eps_site, 0, sd_site, TRUE)) # Pr(random effects)
  lam_i <- rep(0, length(yobs))

  for(i in 1:length(yobs)){
    lam_i[i] = exp(b0 + b1*x1[i] + # fixed effects
                   eps_site[site[i]] # random effects 
                   )
  }
  jnll = jnll - sum(dpois(yobs, lam_i, TRUE)) # Pr(observations|random effects)
  jnll
}

obj = MakeADFun(f, pars, random = "eps_site")
opt = nlminb(obj$par, obj$fn, obj$gr)

sdr = sdreport(obj)
sdr

plot(opt$par[1] + sdr$par.random)
