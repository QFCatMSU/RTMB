# set up unstructured correlation structure for jim...

library(RTMB)

# first, look at unstructured()

unstructured

# now do a little test 

data =  list(k = 3)

pars =  list(
  sigmas = rep(0.1, 3)
)

getAll(data, pars)
us <- unstructured(k) # set up an unstructured object of dim k
us 

OMEGA <- us$corr(sigmas) # correlation matrix 

# Now must set up covariance matrix as 
SIGMA <- diag(sigmas) %*% OMEGA %*% diag(sigmas)

# where SIGMA can be used in dmvnorm
# note sigmas are correlations so must be bounded between correctly

