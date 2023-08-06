library(RTMB)
library(tidyverse)

n_site <- 50
set.seed(12345)
spatial_var <- 1
range <- 0.2
sim1 <- geoR::grf(n = n_site, cov.pars = c(spatial_var, range)) # var, range

x <- sim1$data

n_per_site <- 20

locs <- sim1$coords
n_site <- nrow(locs)

nobs <- n_site * n_per_site
y_obs <- rep(NA, nobs)

beta0 <- 3

i <- 1:n_per_site
for (s in 1:n_site) {
  y_obs[i] <- rpois(1, exp(beta0 + x[s]))
  i <- i + n_per_site
}

data <- list(
  y_obs = y_obs,
  loc_idx = rep(1:n_site, each = n_per_site)
)

mesh <- INLA::inla.mesh.create(locs, refine = TRUE, extend = -0.5, cutoff = 0.01)
spde <- INLA::inla.spde2.matern(mesh, alpha = 2)
data$spde <- spde$param.inla[c("M0", "M1", "M2")]
data$meshidxloc <- mesh$idx$loc # RTMB: 1-based indexing !
n_s <- nrow(data$spde$M0) # Number of points in mesh (including supporting points)

my_df <- data.frame(
  eps = x,
  x = locs[, 1],
  y = locs[, 2]
)

p1 <- my_df %>%
  ggplot(aes(x, y, fill = eps)) +
  geom_point(pch = 21) +
  scale_fill_gradient2() +
  ggtitle("true effects")

parameters <- list(
  beta0 = log(beta0),
  log_tau = -2.0,
  log_kappa = 2.5,
  eps_s = rep(0.0, n_s)
)

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
  Q <- Q_spde(spde, kappa)
  jnll <- 0
  jnll <- jnll - dgmrf(eps_s, 0, Q, log = TRUE, scale = 1 / tau)
  for (i in 1:length(y_obs)) {
    jnll <- jnll - dpois(y_obs[i], exp(beta0 + eps_s[loc_idx[i]]), TRUE)
  }
  range <- sqrt(8) / exp(log_kappa)
  sig_o <- 1 / sqrt(4 * pi * exp(2 * log_tau) * exp(2 * log_kappa))
  ADREPORT(kappa)
  ADREPORT(range)
  ADREPORT(sig_o)
  jnll
}

obj <- MakeADFun(f, parameters, random = "eps_s")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdrep <- sdreport(obj)
truth <- c(sqrt(8)/range, range, sqrt(spatial_var)) 
res <- cbind(sdrep$value, sdrep$value + sdrep$sd %o% c(-1.96, 1.96), truth)
colnames(res) <- c("MLE", "lower95", "upper95", "truth")
print(res)

my_df <- data.frame(
  eps = sdrep$par.random[1:n_site],
  x = locs[, 1],
  y = locs[, 2]
)

p2 <- my_df %>%
  ggplot(aes(x, y, fill = eps)) +
  geom_point(pch = 21) +
  scale_fill_gradient2() +
  ggtitle("estimated random effects")

cowplot::plot_grid(p1, p2)
