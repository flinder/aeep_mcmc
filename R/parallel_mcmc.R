#' Fit Bayesian Models in Parallel
#' 
#' @importFrom foreach foreach %dopar% %do%
#' 
#' @param data A \code{matrix} or \code{data.frame} containing the complete data
#' set.
#' @param cores The number of cores to be used to fit the model. This implies the 
#' number of partitions created. If not supplied the number of cores is detected 
#' automatically.
#' @param combine Method to be used to combine sub-posteriors to the full posterior
#' @param  fun A function fitting the model to the data subsample, first argument 
#' must be dat. No other arguments allowed (for now)
#' (\code{"parametric"}, \code{"semi-parametric"} (or \code{"non-parametric"}))
#' 
#' @return A list containing the full posterior and the sub-posteriors
#' @export

parallel_mcmc <- function(dat, cores, combine = "parametric", fun) {

  # -----------------------------------------
  # Partition the data
  dat <- as.data.frame(dat)
  n_part <- cores
  n <- nrow(dat)
  assgn <- sample(rep(1:n_part, length.out = n))
  pdat <- split(dat, assgn) 
  attr(pdat, "names") <- NULL
    
  # -----------------------------------------
  # Fit model to each subsample
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl) 
  sub_post <- foreach::foreach(i = pdat) %dopar% { fun(dat = i, M = n_part)}
  parallel::stopCluster(cl)
  
  # --------------------------------------------------------
  # Combine sub posteriors to full posterior Return results
  if (combine == "parametric") full_post <- comb_par(sub_post)
  if (combine == "non-parametric") full_post <- comb_npar(sub_post)
  return(list(full = full_post, subs = sub_post))
}
