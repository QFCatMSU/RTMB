library(RTMB)
data <- read.csv("data/sel_dat.csv")

n_sites <- length(unique(data$net_id))
catches <- as.matrix(data[, -c(1, 2)])
colnames(catches) <- NULL

data <- list(
  mesh = c(51, 57, 64, 70, 76, 83, 89, 95, 102, 108, 114, 121, 127),
  n_site = length(unique(data$net_id)),
  n_obs = length(data[, 1]),
  site = rep(
    seq(
      from = 1,
      to = length(unique(data$net_id)),
      length.out = length(unique(data$net_id))
    ),
    each = length(unique(data$length_class))
  ),
  lens = data[, 2],
  catches = catches
)
data$rel_size <- data$mesh / data$mesh[1]

pars <- list(
  ln_k1 = log(150),
  ln_k2 = log(55),
  k1_dev = rep(0, data$n_site),
  k2_dev = rep(0, data$n_site),
  ln_sd_k1 = log(1),
  ln_sd_k2 = log(1)
)

f <- function(pars) {
  getAll(data, pars)
  jnll <- 0
  jnll <- jnll - sum(dnorm(k1_dev, 0, exp(ln_sd_k1), TRUE))
  jnll <- jnll - sum(dnorm(k2_dev, 0, exp(ln_sd_k2), TRUE))
  k1 <- exp(ln_k1 + k1_dev)
  k2 <- exp(ln_k2 + k2_dev)
  sel_mat <- phi_mat <- matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] <- exp(-(lens[i] - k1[site[i]] * rel_size[j])^2 /
        (2 * k2[site[i]]^2 * rel_size[j]^2))
    }
  }
  sel_sums <- rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] <- sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll <- jnll - sum(catches * log(phi_mat))
  jnll
}

obj <- MakeADFun(f, pars, random = c("k1_dev", "k2_dev"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$objective # should solve to 13070.43
sdr <- sdreport(obj)
sdr
