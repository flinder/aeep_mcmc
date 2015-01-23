set_cppo(mode = "fast")
# Example 2pl model

# Simulate some data
Y <- sim_irt(100, 10, seed = 234234)

# Fitting Function
fit <- function(dat, M) {
  Y <- dat
  J <- nrow(Y)
  K <- ncol(Y)
  Y.vec <- stack(as.data.frame(Y))[, 1]
  data.stan <- list("J" = J, 
                    "K" = K, 
                    "N" = J * K,  
                    "jj" = rep(c(1:J), K), 
                    "kk" = rep(c(1:K), each = J), 
                    "Y" = Y.vec,
                    "M" = M
  )
  fileName <- "tests/stan_2pl_dwp.cpp"
  model <- readChar(fileName, file.info(fileName)$size)
  
  res <- rstan::stan(model_code = model, model_name = "2pl", data = data.stan,
                     iter = 5000, warmup = 1000, chains = 1, verbose = FALSE)
  out <- do.call(cbind,res@sim$samples[[1]][- (J + 2 * K + 1)])
  
  # Only do item parameters for now
  out <- out[, -grep("theta", colnames(out))]
  return(out)
}

# Normal estimation
stan_out <- fit(Y, 4)

# Parallelized Estimation
test_par <- parallel_mcmc(dat = Y, cores = 4, combine = "parametric", fun = fit)
