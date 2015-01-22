# Example Linear Regression

# Simulate some data
N <- 1e3
x1 <- rnorm(N, 1, 2)
x2 <- rnorm(N, -1, 2)
y <- x1 - x2 + rnorm(N)
dat <- data.frame(y, x1, x2) 

# Example function: linear regression
fit <- function(data){
  require(rjags)
  N <- nrow(data)
  dataList <- list(x = cbind(data$x1, data$x2), y = data$y, N = N, K = 2)
  initsList <- list(b0 = -1, b = c(0, 0), tau = 1)
  params = c("b0", "b", "tau")
  jmod <- jags.model("linear_regression.txt", data = dataList, n.chains = 1, n.adapt = 500)
  mcmcres <- coda.samples(model = jmod, variable.names = params, 
                          n.iter = 1000, thin = 1)
  return(mcmcres[[1]])
}

out <- parallel_mcmc(data = dat, cores = 4, combine = "parametric")
