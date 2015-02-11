#' Parametric Combination of Posteriors
#' 
#' @param post_list A list of matrices containing samples from the 
#' sub-posteriors
#' @return A list containing a vector of means and the variance-covariance matrix
#' of the full posterior
combine_p <- function(post_list) {
  d <- ncol(post_list[[1]])
  M <- length(post_list)
  var_c <- post_vcm(post_list, d, M)
  mean_c <- post_mean(post_list, var_c, d, M)
  return(list("post_means" = mean_c, "post_variance" = var_c))
}