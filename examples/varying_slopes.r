library(RTMB)
library(tidyverse)
library(MASS)
library(ellipse)
library(ggqfc)
library(gghighlight)

# unstandardized means and sds
n_lakes <- 30 # number of lakes
n_visits <- 10 # number of measurements/years at each lake
total_obs <- n_lakes * n_visits

sigma <- 1 # population (likelihood) sd
betas <- c(14, -0.15) # average intercept and slope
sigmas <- c(2, 1) # intercept and slope sds
rho <- -0.3 # correlation between intercepts and slopes
OMEGA <- matrix(c(1, rho, rho, 1), nrow = 2) # correlation matrix
SIGMA <- diag(sigmas) %*% OMEGA %*% diag(sigmas) # var covar

# draw correlated slopes and intercepts:
set.seed(5621)
vary_effects <- mvrnorm(n_lakes, betas, SIGMA)

b0 <- vary_effects[, 1]
b1 <- vary_effects[, 2]

# visualize it
plot(b0, b1,
  col = "steelblue", xlab = "intercepts (b0)", ylab = "slopes(b1)",
  ylim = c(-3, 3), xlim = c(10, 19), main = "Visualizing correlation between slopes and intercepts"
)
for (l in c(0.1, 0.3, 0.5, 0.8, 0.99)) {
  lines(ellipse(SIGMA, centre = betas, level = l))
}

# simulating the rest of our data
n <- n_visits * n_lakes # total number of observations
visit <- rep(1:n_visits, n_lakes)
lake_id <- rep(1:n_lakes, each = n_visits) # create a lake ID
x <- rnorm(n, 0, 1) # create fake covariate data
mu <- b0[lake_id] + b1[lake_id] * x # mean prediction
y <- rnorm(n, mu, sigma) # add likelihood error to mean prediction

sim_data <- data.frame(y, x, lake_id, visit)

# say bad weather means you couldn't go out 25% of the time:
keep <- rbinom(nrow(sim_data), size = 1, prob = 0.75)
sim_data <- sim_data[which(keep == 1), ]

p1 <-
  sim_data %>%
  ggplot(aes(x = x, y = y, group = lake_id)) +
  geom_point(color = "firebrick", alpha = 0.50) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  xlab("standardized juvenile density") +
  ylab("growth rate cm/yr") +
  facet_wrap(~lake_id) +
  theme_qfc()
p1

p2 <- sim_data %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(color = "firebrick", alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  xlab("standardized juvenile density") +
  ylab("growth rate cm/yr") +
  theme_qfc()
p2

p3 <-
  sim_data %>%
  ggplot(aes(x = x, y = y)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_point(size = 3, color = "firebrick", alpha = 0.5) +
  gghighlight(lake_id == 18,
    use_group_by = FALSE, max_highlight = Inf,
    use_direct_label = FALSE
  ) +
  labs(
    subtitle = "Red points are data from lake 18 \n
       Grey points are data from all other lakes"
  ) +
  xlab("juvenile density") +
  ylab("growth rate cm/yr") +
  theme_qfc()
p3

#--------------------------------
# estimate it using RTMB

data <- list(
  y = sim_data$y,
  x = sim_data$x,
  lake = sim_data$lake_id, 
  nobs = nrow(sim_data)
)

pars <- list(
  global_betas = c(13, 0), # overall average estiamtes
  trans_rho = -0.4,
  log_sd_obs = 0.2,
  log_sds_re = c(-1, -1),
  betas_lake = matrix(0, ncol = 2, nrow = n_lakes) # lake specific deviations
)

to_cor <- function(x) {
  2 / (1 + exp(-2 * x)) - 1
}

f <- function(pars) {
  getAll(data, pars)
  rho <- to_cor(trans_rho)
  sd_obs <- exp(log_sd_obs)
  sds_re <- exp(log_sds_re)
  jnll <- 0
  OMEGA <- matrix(c(1, rho, rho, 1), nrow = 2) # correlation matrix
  SIGMA <- diag(sds_re) %*% OMEGA %*% diag(sds_re) # calculate var covar from correlation matrix 
  jnll <- jnll - sum(dmvnorm(betas_lake, global_betas, SIGMA, TRUE)) # MVN re's
  y_pred <- rep(0, nobs)
  for (i in 1:nobs) {
    y_pred[i] <- betas_lake[lake[i], 1] + betas_lake[lake[i], 2] * x[i] # mean prediction
  }
  jnll = jnll - sum(dnorm(y, y_pred, sd_obs, TRUE))
  REPORT(rho)
  jnll
}

obj <- MakeADFun(f, pars, random = "betas_lake")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr = sdreport(obj)

sdr
obj$report()

# pluck the random effects out of the sdreport
res <- sdr$par.random

# vs. truth
truth <- c(b0, b1)

dim(truth) <- dim(res) <- c(n_lakes, 2)

truth
res
