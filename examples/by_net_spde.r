library(RTMB)
library(tidyverse)

data = read.csv("data/sel_dat.csv")
n_sites = length(unique(data$net_id))
catches = as.matrix(data[, which(grepl("mesh", colnames(data)))])
colnames(catches) = NULL

locs <- unique(data[, c("easting_km", "northing_km")])
plot(locs)

mesh = INLA::inla.mesh.2d(loc=locs, max.edge=c(6,50))

plot(mesh)
points(locs, col = "#075057", pch = 16)
mesh$n
n_sites

data <- within(data, net <- as.numeric(interaction(data$net_id, drop = TRUE, lex.order = F)))
data <- data[order(data$net), ]

data <- list(
  meshes = c(51, 57, 64, 70, 76, 83, 89, 95, 102, 108, 114, 121, 127),
  n_site = length(unique(data$net)),
  n_obs = nrow(data),
  loc_idx = data$net, 
  lens = data[, 2],
  catches = catches
)
data$rel_size = data$meshes / data$meshes[1]

spde <- INLA::inla.spde2.matern(mesh, alpha = 2)

data$spde <- spde$param.inla[c("M0", "M1", "M2")]

pars = list(
  ln_k1 = 5.6,
  ln_k2 = 4.2,
  eps_k1 = rep(0, mesh$n),
  eps_k2 = rep(0, mesh$n),
  log_tau1 = 2.12,
  log_tau2 = -2.0, 
  log_kappa1 = -1.43, 
  log_kappa2 = 0
)

Q_spde <- function(spde, kappa) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2
}

spde <- INLA::inla.spde2.matern(mesh, alpha = 2)
spdeMatrices <- spde$param.inla[c("M0", "M1", "M2")]

f <- function(parameters) {
  getAll(data, parameters)
  tau1 <- exp(log_tau1)
  tau2 <- exp(log_tau2)
  kappa1 <- exp(log_kappa1) 
  kappa2 <- exp(log_kappa2) 
  
  Q1 <- Q_spde(spde, kappa1)
  Q2 <- Q_spde(spde, kappa2)

  jnll <- 0
  jnll <- jnll - dgmrf(eps_k1, 0, Q1, log = TRUE, scale = 1 / tau1) # mean0 k1 deviates
  jnll <- jnll - dgmrf(eps_k2, 0, Q2, log = TRUE, scale = 1 / tau2) # mean0 k2 deviates
  
  k1 = exp(ln_k1 + eps_k1)
  k2 = exp(ln_k2 + eps_k2)
  
  sel_mat = phi_mat = matrix(0, nrow(catches), ncol(catches))

  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] = exp(-(lens[i] - k1[loc_idx[i]] * rel_size[j])^2 /
        (2 * k2[loc_idx[i]]^2 * rel_size[j]^2))
    }
  }
  sel_sums = rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] = sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll = jnll - sum(catches * log(phi_mat))

  # Derived variables 
  range1 <- sqrt(8) / exp(log_kappa1)
  range2 <- sqrt(8) / exp(log_kappa2)

  sig_o_k1 <- 1 / sqrt(4 * pi * exp(2 * log_tau1) * exp(2 * log_kappa1))
  sig_o_k2 <- 1 / sqrt(4 * pi * exp(2 * log_tau2) * exp(2 * log_kappa2))

  ADREPORT(range1)
  ADREPORT(range2)
  ADREPORT(sig_o_k1)
  ADREPORT(sig_o_k2)

  jnll
}

obj <- MakeADFun(f, pars, random = c("eps_k1", "eps_k2"))
obj$fn()
obj$gr()

opt <- nlminb(obj$par, obj$fn, obj$gr)
sdrep <- sdreport(obj)
sdrep
