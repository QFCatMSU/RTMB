library(RTMB)
library(mvtnorm)

#-------------------------------------------------------------------------------------
# A spatial Poisson model, loosely following coding style in  
# Anderson and Ward 2019: Black swans in space
# NOTE: this model currently does NOT account for spatiotemporal extremes

# model form:
# y_i ~ normal(mu_i)
# mu_i = exp(b0 + eps(site(i)) )
# eps(s) ~ GRF(0, SIGMA)
# SIGMA  = gp_sigma^2 * exp(-dist_sites/gp_theta)
#
# where gp_sigma, gp_theta are estimated Gaussian process parameters
# and dist_sites is a distance matrix indicating euclidean distances among sites

# Note that we need to simulate sites, years, and a distance matrix

#-------------------------------------------------------------------------------------
# simulate it 
set.seed(17)
n_site <- 50 # number of sampling locations for each year
n_per_site <- 30
# simulate random x,y site locations:
g <- data.frame(
  lon = runif(n_site, 0, 10),
  lat = runif(n_site, 0, 10)
)

locs <- unique(g)
dist_sites <- as.matrix(dist(locs))

# set up a vector indicating the site index identifiers
site_id <- rep(1:n_site, each = n_per_site)
nobs = length(site_id)
# model parameters to simulate
gp_theta <- 1 # Gaussian process scale parameter
gp_sigma <- 0.15 # Gaussian process variance / spatial noise parameter
beta0 <- log(1.75) # linear regression intercept

gp_sigma_sq = gp_sigma*gp_sigma
SIGMA = gp_sigma_sq * exp(- dist_sites / gp_theta)

eps_s = rmvnorm(nobs, rep(0, n_site), SIGMA)

