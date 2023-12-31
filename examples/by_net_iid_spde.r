library(RTMB)
library(tidyverse)

data <- read.csv("data/sel_dat.csv")
year <- data$year - min(data$year) + 1

n_year <- 7 # 2007 to 2013

n_sites <- length(unique(data$net_id))
catches <- as.matrix(data[, which(grepl("mesh", colnames(data)))])
colnames(catches) <- NULL

locs <- unique(data[, c("easting_km", "northing_km")])

mesh <- INLA::inla.mesh.2d(loc = locs, max.edge = c(6, 50))

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
  year_idx = year,
  n_year = n_year,
  lens = data[, 2],
  catches = catches
)
data$rel_size <- data$meshes / data$meshes[1]

spde <- INLA::inla.spde2.matern(mesh, alpha = 2)

data$spde <- spde$param.inla[c("M0", "M1", "M2")]

pars <- list(
  ln_k1 = 6.69,
  ln_k2 = 4.5,
  eps_k1_st = matrix(0.0, nrow = mesh$n, ncol = n_year),
  eps_k2_st = matrix(0.0, nrow = mesh$n, ncol = n_year),
  log_tau1 = 1.02,
  log_tau2 = 1.02,
  log_kappa1 = -.4,
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
  catches <- OBS(catches)
  tau1 <- exp(log_tau1)
  tau2 <- exp(log_tau2)
  kappa1 <- exp(log_kappa1)
  kappa2 <- exp(log_kappa2)
  Q1 <- Q_spde(spde, kappa1)
  Q2 <- Q_spde(spde, kappa2)

  jnll <- 0
  for (t in 1:n_year) {
    jnll <- jnll - dgmrf(eps_k1_st[, t], 0, Q1, TRUE, scale = 1 / tau1)
    jnll <- jnll - dgmrf(eps_k2_st[, t], 0, Q2, TRUE, scale = 1 / tau2)
  }

  k1 <- exp(ln_k1 + eps_k1_st)
  k2 <- exp(ln_k2 + eps_k2_st)

  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))

  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1[loc_idx[i], year_idx[i]] * rel_size[j])^2 /
        (2 * k2[loc_idx[i], year_idx[i]]^2 * rel_size[j]^2))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  # jnll = jnll - sum(catches * log(phi_mat))
  jnll <- jnll - sum(dpois(catches, phi_mat, log = TRUE)) # coding to allow simulation

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

obj <- MakeADFun(f, pars, random = c("eps_k1_st", "eps_k2_st"))
obj$fn()
obj$gr()

opt <- nlminb(obj$par, obj$fn, obj$gr)
sdrep <- sdreport(obj)
sdrep

sdrep$value
sdrep$sd

# set.seed(6266)
# system.time(
#   chk <- checkConsistency(obj, estimate = TRUE, n = 500)
# )

# summary(chk)$estimate$par
# par(mfrow = c(2, 3))
# hist(summary(chk)$estimate$par$ln_k1, breaks = 20)
# abline(v = summary(chk)$par[1])

# hist(summary(chk)$estimate$par$ln_k2, breaks = 20)
# abline(v = summary(chk)$par[2])

# hist(summary(chk)$estimate$par$log_tau1, breaks = 20)
# abline(v = summary(chk)$par[3])

# hist(summary(chk)$estimate$par$log_tau2, breaks = 20)
# abline(v = summary(chk)$par[4])

# hist(summary(chk)$estimate$par$log_kappa1, breaks = 20)
# abline(v = summary(chk)$par[5])

# hist(summary(chk)$estimate$par$log_kappa2, breaks = 20)
# abline(v = summary(chk)$par[6])
