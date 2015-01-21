# ======================================================================
# Function for parallel sampling from bayesian posteriors as described in 
# Neiswanger et al. (2013)
# 
# Fridolin Linder
#
# Arguments:
# data: original data set
# cores: number of cores to be used = number of sub posteriors estimated
# combine: Parametric, Semiparametric or Non- Parametric combination of
#          sub posteriors (Semi not implemented yet and non very slow)
# fit: A function fitting the model to the data subsample, first argument 
#      must be data. No other arguments allowed (for now)

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
  cl <- makeCluster(cores)
  registerDoParallel(cl) 
  
  # Fit models
  sub_post <- foreach(i = pdat, .packages = "MCMCpack") %dopar% {
    fit(data = i)
  }
  stopCluster(cl)
  
  # -----------------------------------------
  # Combine sub posteriors to full posterior
  # ----------------------------------------
  
  # Parametric combination of subposteriors
  # mcmcout: list of matrices containing posterior samples from the subposteriors
  comb_par <- function(mcmcout) {
    vcms <- lapply(mcmcout, var)
    ivcms <- lapply(vcms, solve)
    var.c <- solve(Reduce("+", ivcms))
    means <- lapply(mcmcout, function(x) apply(x, 2, mean))
    mprod <- mapply(function(x, y) x %*% y, ivcms, means, SIMPLIFY = F)
    mean.c <- var.c %*% Reduce("+", mprod)
    return(list("post_means" = mean.c, "post_variance" = var.c))
  }
  
  # Non-parametric combination of subposteriors
  # mcmcout: list of matrices containing posterior samples from the subposteriors
  comb_npar <- function(mcmcout) {
    require(mvtnorm)
    d <- ncol(mcmcout[[1]])
    T_ <- nrow(mcmcout[[1]])
    M <- length(mcmcout)
    
    w <- function(ind, mcmcout, h, d, mn = FALSE) {
      ind <- as.list(ind)
      sel <- mapply(function(x, y) x[y, ], mcmcout, ind)
      theta_bar <- apply(sel, 1, mean)
      if(mn) return(theta_bar)
      out <- prod(apply(sel, 2, function(x) dmvnorm(x, theta_bar, h^2 * diag(d))))
      return(out)
    }
    
    t_dot <- sample(c(1:T_), M, replace = T)
    c_dot <- t_dot
    
    out <- matrix(NA, ncol = d, nrow = T_)  
    for(i in 1:T_) {
      h <- i^(- 1 / (4 + d))
      for(m in 1:M) {
        c_dot <- t_dot
        c_dot[m] <- sample(c(1:T_), 1)
        if(runif(1) < w(c_dot, mcmcout, h, d) / w(t_dot, mcmcout, h, d)) {
          t_dot <- c_dot
        } 
      }
      out[i, ] <- rmvnorm(1, w(t_dot, mcmcout, h, d, mn = T), h^2/M * diag(d))
    }
    return(out)
  }
  
  # -----------------------------------------
  # Return results
  # ----------------------------------------
  
  if (combine == "parametric") full_post <- comb_par(sub_post)
  if (combine == "non parametric") full_post <- comb_npar(sub_post)
  return(list(full = full_post, subs = sub_post))
}
