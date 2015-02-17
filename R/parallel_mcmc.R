#' Fit Bayesian Models in Parallel
#' 
#' @importFrom foreach foreach %dopar% %do%
#' 
#' @param data A matrix or data.frame containing the complete data
#' set.
#' @param cores The number of cores to be used to fit the model. This implies the 
#' number of partitions created. If not supplied the number of cores is detected 
#' automatically.
#' @param combine Method to be used to combine sub-posteriors to the full posterior
#' (\code{"parametric"}, \code{"semi-parametric"} (or \code{"non-parametric"}))
#' @param  fun A function fitting the model to the data subsample. For now the
#' function must only have two arguments. The first beeing the data, the second 
#' beeing the number of partitions.
#' 
#' @return A list containing the full posterior and the sub-posteriors
#' @export

parallel_mcmc <- function(data, cores = NULL, combine = NULL, fun, ...) {

  # -----------------------------------------
  # Partition the data
  data <- as.data.frame(data)
  if(is.null(cores)) cores <- detectCores()
  n_part <- cores
  n <- nrow(data)
  assgn <- sample(rep(1:n_part, length.out = n))
  pdat <- split(data, assgn) 
  attr(pdat, "names") <- NULL
  
  # -----------------------------------------
  # Fit model to each subsample
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl) 
  sub_post <- foreach::foreach(i = pdat, .packages = "rstan") %dopar% { 
    fun(data = i, n_partitions = n_part, ...)
    #fun(data = i, n_partitions = n_part, c_mod = c_model_dw, iter = iter, chains = chains, 
    #    warmup = warmup)
  }
  parallel::stopCluster(cl)
  
  # --------------------------------------------------------
  # Combine sub posteriors to full posterior Return results
  if(is.null(combine)) {
    out <- sub_post
  } else if(combine == "parametric") {
    out <- combine_p(sub_post)
  } else if(combine == "non-parametric") {
    out <- combine_np(sub_post)  
  } else if(combine == "semi-parametric") {
    out <- combine_sp(sub_post)  
  } else {
    warning("Invalid 'combine' argument, returning sub-posteriors.
            Sub-posteriors can be combined with combine_p(), combine_sp() and
            combine_np()")
    out <- sub_post
  }
  return(out)
}
