#' Fit 2pl IRT Models in Parallel
#' 
#' @param data A matrix or data.frame containing the complete data
#' set.
#' @param cores The number of cores to be used to fit the model. This implies the 
#' number of partitions created. If not supplied the number of cores is detected 
#' automatically.
#' @param iter Number of iterations for each mcmc sampler.
#' @param warmup Number of warmup samples for each sampler.
#' @param chains Number of chains to run on each sub-sample.
#' @param combine Method to be used to combine sub-posteriors to the full posterior
#' (\code{"parametric"}, \code{"semi-parametric"} (or \code{"non-parametric"}))
#' 
#' @return A list containing the full posterior and the sub-posteriors.
#' 
#' If \code{combine == "parametric"} The vector of means and the variance-covariance
#' matrix of the full-posterior is returned
#' 
#' If \code{combine == "non-parametric"} The full posterior is returned as a 
#' matrix with number of rows equal to number of iterations minus warmup and 
#' number of columns equal to the number of free parameters in the model.
#' 
#' Note: For now the function only returns the item parameters. INdividual level
#' parameters will be implemented shortly.
#' @export

par2pl <- function(data, 
                   cores = NULL, 
                   iter, 
                   warmup = floor(iter/2), 
                   chains,
                   combine = "parametric") {
  
  if(is.null(cores)) cores = parallel::detectCores()
  
  # Fitting function
  fit <- function(data, M) {
    J <- nrow(data)
    K <- ncol(data)
    data <- as.data.frame(data)
    class(data) <- "list"
    data_stan <- append(data, list(J = J, K = K, M = M))

    # Create stan model string
    
    model <- paste0("
    functions {
      real dw_normal_log(real y, real mu, real sigma, int M) {
        return normal_log(y, mu, sigma) / M;
      }
      real dw_lognormal_log(real y, real mu, real sigma, int M) {
        return lognormal_log(y, mu, sigma) / M;
      }
    }
    data {
      int<lower=1> J;
      int<lower=1> K;
      int<lower=1> M;",
    paste0("int<lower=0, upper=1> V", c(1:K), "[J];", collapse = "\n "),
    "}
    parameters {
      real<lower=0> alpha[K]; 
      real beta[K]; 
      vector[J] theta; 
    }
    model {
      for(j in 1:J) {
        theta[j] ~ dw_normal(0.0, 1.0, M);
      }
      for(k in 1:K) {
        beta[k] ~ dw_normal(0.0, 1.0, M);
        alpha[k] ~ dw_lognormal(0.0, 6, M);
      }",
    paste0("V", c(1:K), " ~ bernoulli_logit(beta[", c(1:K), "] + alpha[", c(1:K)
           , "] * theta);", collapse = "\n "),
    "}
    ")
   
    res <- rstan::stan(model_code = model, data = data_stan, iter = iter, 
                       warmup = warmup, chains = chains)
    out <- do.call(cbind,res@sim$samples[[1]][- (J + 2 * K + 1)])
    
    # Only do item parameters for now
    out <- out[, -grep("theta", colnames(out))]
    return(out)
  }
  
  # Fit in parallel
  out <- parallel_mcmc(data = Y, cores = cores, combine = combine, fun = fit)
  
  return(out)
}
