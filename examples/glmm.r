library(RTMB)

#-------------------------------------------------------------------
# simulate and estimate a Poisson GLMM
# y_i ~ Poisson(lambda_i)
# log(lambda_i) = b0 + b1 + eps_site
# where eps_site ~ N(0, sd_eps)

#--------------------------------------------------------------------

b0 = 3
b1 = 2.3
sd_eps = 0.7
n_site = 15
n_replicates = 5
site = rep(1:n_replicates, each = n_site)
nobs = length(site)

set.seed(333)
x1 = runif(nobs, -1, 1)
eps_i = rnorm(n_site, 0, sd_eps)
lambda_i = rep(NA, n_site * n_replicates)
nobs = length(lambda_i)
lambda_i = exp(b0 + b1 * x1 + eps_i)
yobs = rpois(nobs, lambda_i)

#-------------------------------------------------------------------
# estimate it

data = list(
  yobs = yobs,
  Xmat = model.matrix(~ 1 + x1)
)

pars = list(
  bvec = c(1, 1), # b0, b1
  log_sd_site = log(1),
  eps_site = rep(0, length(data$yobs))
)

f = function(pars) {
  getAll(data, pars)
  sd_site = exp(log_sd_site)                           # back transform
  jnll = 0                                             # initialize
  jnll = jnll - sum(dnorm(eps_site, 0, sd_site, TRUE)) # Pr(random effects)
  lam_i = exp(                                         # link f(x)
    Xmat %*% bvec +                                    # fixed effects
    eps_site                                           # random effects
  )
  jnll = jnll - sum(dpois(yobs, lam_i, TRUE))          # Pr(observations|random effects)
  jnll
}

obj = MakeADFun(f, pars, random = "eps_site")
opt = nlminb(obj$par, obj$fn, obj$gr)

sdr = sdreport(obj)
sdr

