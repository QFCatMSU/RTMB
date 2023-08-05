library(RTMB)
library(mvtnorm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NOTE-this is not yet correct and needs to be fixed...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#--------------------------------------------------------------------------------
# A hierarchical Poisson model with a spatially-explicit random effect, 
# loosely following coding style/math in Anderson and Ward 2019 Ecology

# model form:
# y_i ~ pois(lambda_i)
# mu_i = exp(b0 + eps(site(i)) )
# eps(s) ~ GRF(0, SIGMA)
# SIGMA  = gp_sigma^2 * exp(-dist_sites/gp_theta)
#
# where gp_sigma, gp_theta are estimated Gaussian process parameters
# and dist_sites is a distance matrix indicating euclidean distances among sites

# Note that we need to simulate sites, years, and a distance matrix

#--------------------------------------------------------------------------------
# simulate it

set.seed(17)
n_site = 50 # number of sampling locations for each year
n_per_site = 30
# simulate random x,y site locations:
g = data.frame(
  lon = runif(n_site, 0, 10),
  lat = runif(n_site, 0, 10)
)

locs = unique(g)
dist_sites = as.matrix(dist(locs))

# set up a vector indicating the site index identifiers
site_id = rep(1:n_site, each = n_per_site)
nobs = length(site_id)

# model parameters to simulate
gp_theta = 1      # Gaussian process scale parameter
gp_sigma = 0.15   # Gaussian process variance / spatial noise parameter
beta0 = log(1.75) # linear regression intercept

gp_sigma_sq = gp_sigma * gp_sigma
SIGMA = gp_sigma_sq * exp(-dist_sites / gp_theta)

set.seed(666)
eps_s = rmvnorm(nobs, rep(0, n_site), SIGMA)
y_det = y_obs = rep(0, nobs)
y_det = exp(beta0 + eps_s)
y_obs = rpois(nobs, y_det)

#-------------------------------------------------------------------------------------
# estimate it 

data = list(
  y_obs = y_obs,
  dist_sites = dist_sites
)

pars = list(
  beta0 = log(1.75),
  gp_theta = 0.7, 
  gp_sigma = 0.1,
  eps_s = rep(0, n_site)
)

f = function(pars){
  getAll(data, pars)
  SIGMA = gp_sigma * exp(-dist_sites / gp_theta)
  jnll = 0
  jnll = jnll - sum(dmvnorm(eps_s, SIGMA, TRUE))
  y_hat = exp(beta0 + eps_s) # index on site_i in more complex examples
  jnll = jnll - sum(dpois(y_obs, y_hat), TRUE)
  jnll 
}

obj = MakeADFun(f, pars, random = "eps_s")
opt = nlminb(obj$par, obj$fn, obj$gr)
opt

sdr = sdreport(obj)
sdr
