library(RTMB)
library(tidyverse)

# generate a spatial-temporal field evolving via AR-1
set.seed(666)
n_site <- 30
n_year <- 10
omega_dev_st <- matrix(0, nrow = n_site, ncol = n_year)
spatial_var <- 1
range <- 0.2
rho <- 0.5

# intialize, get locs
sim <- geoR::grf(n = n_site, cov.pars = c(spatial_var, range)) # var, range
omega_dev_st[, 1] <- sim$data
locs <- sim$coords

for (t in 2:n_year) {
  sim <- geoR::grf(n = n_site, grid = locs, cov.pars = c(spatial_var, range)) # var, range
  omega_dev_st[, t] <- sim$data
}

eps_st <- omega_dev_st
for (t in 2:n_year) {
  eps_st[, t] <- rho * eps_st[, t - 1] + sqrt(1 - rho^2) * omega_dev_st[, t]
}

my_df <- data.frame(
  eps = c(eps_st),
  x = rep(locs[, 1], n_year),
  y = rep(locs[, 2], n_year),
  year = rep(1:n_year, each = n_site),
  site = rep(1:n_site, n_year)
)

p1 <- my_df %>%
  ggplot(aes(x, y, fill = eps)) +
  geom_point(pch = 21) +
  scale_fill_gradient2() +
  facet_wrap(~year) +
  ggtitle(expression(true ~ epsilon[st]))
p1

# generate observations based on eps_st
n_per_site <- 20
n_obs <- n_year * n_site * n_per_site
y_obs <- rep(NA, n_obs)

beta0 <- 3
i <- 1:n_per_site
for (s in 1:n_site) {
  for (y in 1:n_year) {
    y_obs[i] <- rpois(n_per_site, exp(beta0 + eps_st[s, t]))
    i <- i + n_per_site
  }
}

data <- list(
  y_obs = y_obs,
  site = rep(1:n_site, each = n_per_site * n_year),
  year = rep(1:n_year, each = n_per_site * n_site)
)

mesh <- INLA::inla.mesh.create(locs,
  refine = TRUE,
  extend = -0.5, cutoff = 0.01
)
spde <- INLA::inla.spde2.matern(mesh, alpha = 2)
data$spde <- spde$param.inla[c("M0", "M1", "M2")]
n_s <- nrow(data$spde$M0) # Number of points in mesh (including supporting points)

parameters <- list(
  beta0 = log(beta0),
  log_tau = -2.0,
  log_kappa = 2.5,
  trho = 0.5,
  eps_st = matrix(0.0, nrow = n_s, ncol = n_year)
)

to_cor <- function(x) {
  2 / (1 + exp(-2 * x)) - 1
}

# Objective function
Q_spde <- function(spde, kappa) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2
}

f <- function(parameters) {
  getAll(data, parameters)
  tau <- exp(log_tau)
  kappa <- exp(log_kappa)
  rho <- to_cor(trho)
  Q <- Q_spde(spde, kappa)
  jnll <- 0
  jnll <- jnll - dgmrf(eps_st[, 1], 0, Q, TRUE) # initialize
  for (t in 2:n_year) { # autoregress
    jnll <- jnll - dgmrf(eps_st[, t], rho * eps_st[, t - 1], Q, TRUE, sqrt(1 - rho^2))
  }
  for (i in 1:length(y_obs)) {
    jnll <- jnll - dpois(y_obs[i], exp(beta0 + eps_st[site[i], year[i]] / tau), TRUE)
  }
  range <- sqrt(8) / exp(log_kappa)
  sig_o <- 1 / sqrt(4 * pi * exp(2 * log_tau) * exp(2 * log_kappa))
  ADREPORT(kappa)
  ADREPORT(range)
  ADREPORT(sig_o)
  ADREPORT(rho)
  jnll
}

obj <- MakeADFun(f, parameters, random = "eps_st")
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$convergence
sdrep <- sdreport(obj)
truth <- c(sqrt(8) / range, range, sqrt(spatial_var), rho)
res <- cbind(sdrep$value, sdrep$value + sdrep$sd %o% c(-1.96, 1.96), truth)
colnames(res) <- c("MLE", "lower95", "upper95", "truth")
print(res)

my_df <- data.frame(
  eps = c(sdrep$par.random[1:n_site]),
  x = rep(locs[, 1], n_year),
  y = rep(locs[, 2], n_year),
  year = rep(1:n_year, each = n_site),
  site = rep(1:n_site, n_year)
)

p2 <- my_df %>%
  ggplot(aes(x, y, fill = eps)) +
  geom_point(pch = 21) +
  scale_fill_gradient2() +
  facet_wrap(~year) +
  ggtitle(expression(estimated ~ epsilon[st]))
p2

# cowplot::plot_grid(p1, p2)
