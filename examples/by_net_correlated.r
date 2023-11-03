library(RTMB)
data = read.csv("data/sel_dat.csv")

n_sites = length(unique(data$net_id))
catches = as.matrix(data[, which(grepl("mesh", colnames(data)))])
colnames(catches) = NULL

data = list(
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
data$rel_size = data$mesh / data$mesh[1]

pars = list(
  log_global_k = c(log(150), log(55)),
  k_site = matrix(0, ncol = 2, nrow = data$n_site), # matrix of k1, k2 values, with nrow = n_site and col 1 = k1 
  log_sds_re = c(log(1), log(1)), 
  trans_rho = 0
)
to_cor <- function(x) {
  2 / (1 + exp(-2 * x)) - 1
}

f = function(pars) {
  getAll(data, pars)
  sds_re = exp(log_sds_re)
  rho = to_cor(trans_rho)
  jnll = 0
  OMEGA <- matrix(c(1, rho, rho, 1), nrow = 2) # correlation matrix
  SIGMA <- diag(sds_re) %*% OMEGA %*% diag(sds_re) # calculate var covar from correlation matrix 
  jnll <- jnll - sum(dmvnorm(k_site, c(0,0), SIGMA, TRUE)) # MVN re's
  k_mat <- matrix(0, ncol = 2, nrow = n_site)
  k_mat[,1] <- exp(log_global_k[1] + k_site[,1]) 
  k_mat[,2] <- exp(log_global_k[2] + k_site[,2])
  sel_mat = phi_mat = matrix(0, nrow(catches), ncol(catches))
  for (i in 1:nrow(sel_mat)) {
    for (j in 1:ncol(sel_mat)) {
      sel_mat[i, j] = exp(-(lens[i] - k_mat[site[i], 1] * rel_size[j])^2 /
        (2 * k_mat[site[i], 2]^2 * rel_size[j]^2))
    }
  }
  sel_sums = rowSums(sel_mat)
  for (i in 1:nrow(phi_mat)) {
    for (j in 1:ncol(phi_mat)) {
      phi_mat[i, j] = sel_mat[i, j] / sel_sums[i]
    }
  }
  jnll = jnll - sum(catches * log(phi_mat))
  jnll
}

obj = MakeADFun(f, pars, random = c("k_site"))
opt = nlminb(obj$par, obj$fn, obj$gr)
opt$objective 
sdr = sdreport(obj)
sdr
