#' Parametric Combination of Posteriors
#' 
#' @param post_list A \code{list} of matrices containing samples from the 
#' sub-posteriors
#' @return A list containing a vector of means and variances of the full 
#' posterior

comb_par <- function(post_list) {
  vcms <- lapply(post_list, var)
  ivcms <- lapply(vcms, solve)
  var.c <- solve(Reduce("+", ivcms))
  means <- lapply(post_list, function(x) apply(x, 2, mean))
  mprod <- mapply(function(x, y) x %*% y, ivcms, means, SIMPLIFY = F)
  mean.c <- var.c %*% Reduce("+", mprod)
  return(list("post_means" = mean.c, "post_variance" = var.c))
}
