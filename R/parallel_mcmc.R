#' Fit Bayesian Models in Parallel
#' 
#' @param data A \code{matrix} or \code{data.frame} containing the complete data
#' set.
#' @param cores The number of cores to be used to fit the model. This implies the 
#' number of partitions created. If not supplied the number of cores is detected 
#' automatically.
#' @param combine Method to be used to combine sub-posteriors to the full posterior
#' @param  fit A function fitting the model to the data subsample, first argument 
#' must be data. No other arguments allowed (for now)
#' (\code{"parametric"}, \code{"semi-parametric"} (or \code{"non-parametric"}))
#' @return A list containing the full posterior and the sub-posteriors

parallel_mcmc <- function(data, cores, combine = "parametric", fit) {

  ## Partition the data
  data <- as.data.frame(data)
  n_part <- cores
  n <- nrow(data)
  assgn <- sample(rep(1:n_part, length.out = n))
  pdat <- split(data, assgn) 
  attr(pdat, "names") <- NULL
    
  # -----------------------------------------
  # Fit model to each subsample
  # ----------------------------------------
  
  # Register Workers for parallel computing
  cl <- doParallel::makeCluster(cores)
  doParallel::registerDoParallel(cl) 
  
  # Fit models
  sub_post <- doParallel::foreach(i = pdat, .packages = "MCMCpack") %dopar% {
    fit(data = i)
  }
  doParallel::stopCluster(cl)
  
  # --------------------------------------------------------
  # Combine sub posteriors to full posterior Return results
  # --------------------------------------------------------
  
  if (combine == "parametric") full_post <- comb_par(sub_post)
  if (combine == "non parametric") full_post <- comb_npar(sub_post)
  return(list(full = full_post, subs = sub_post))
}
