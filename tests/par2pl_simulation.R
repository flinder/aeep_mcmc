## Simulation for testing par2pl.R compared to normal fit and subsample method
library(doParallel)
library(rstan)
library(irtpar)

fit_stan <- function(data, iter, warmup, chains) {
  stan_model <- "
        data {
          int<lower=1> J; // number of students
          int<lower=1> K; // number of questions
          int<lower=1> N; // number of observations
          int<lower=1,upper=J> jj[N]; // student for observation n
          int<lower=1,upper=K> kk[N]; // question for observation n
          int<lower=0,upper=1> Y[N]; // correctness for observation n
        }
        parameters {
          real<lower=0, upper=5> alpha[K]; // discrimination parameters
          real beta[K]; // difficulty parameters
          real theta[J]; // student ability parameters
        }
        model {
          real pi[N];
          alpha ~ lognormal(0, 6); 
          beta ~ normal(0, 1); 
          theta ~ normal(0, 1); 
          for (n in 1:N)
            pi[n] <- 1 / (1 + exp( - alpha[kk[n]] * (theta[jj[n]] - beta[kk[n]])));
    
          Y ~ bernoulli(pi);
        }    
    "
  Y <- data
  J <- nrow(Y)
  K <- ncol(Y)
  Y_vec <- stack(as.data.frame(Y))[, 1]
  data_stan <- list("J" = J, 
                    "K" = K, 
                    "N" = J * K,  
                    "jj" = rep(c(1:J), K), 
                    "kk" = rep(c(1:K), each = J), 
                    "Y" = Y_vec
                    )
  rstan::stan(model_code = stan_model, data = data_stan, iter = iter, 
              warmup = warmup, chains = chains)
}

### Simulation:

# Parameters
MC_iter <- 50
n_subjects <- 1e2
n_items <- 20
seed <- 277083
cores_f <- 4

for(i in 1:MC_iter) {
  # Simulate data
  Y <- sim_irt(n_subjects, n_items, seed = seed)
  
  # Estimate Full model (3 Chains in parallel)
  ptm <- proc.time()
  #cl <- makeCluster(cores_f)
  #registerDoParallel(cl) 
  fitlist <- foreach(i = 1:3, 
                     .combine = sflist2stanfit, 
                     .packages = "rstan", export = "fit_stan") %dopar% {
                       fit_stan(Y, iter = 1000, warmup = 500, chains = 3)
                     }
  #stopCluster(cl)
  proc.time() - ptm
  
  foreach(i = list(1,2,3)) %do% i  
  
  # Estimate sub sample model
  
  
  # Estimate Parallel model
  
}




