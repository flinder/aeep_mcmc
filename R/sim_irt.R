#' Simulate IRT Data
#' 
#' @importFrom foreach foreach %dopar% %do%
#' 
#' @param n_subjects Number of subjects (rows).
#' @param n_items Number of items (columns).
#' @param p_theta A vector of length 2 containing mean and standard deviation for
#' the normal distribution the thetas are generated from.
#' @param p_alpha A vector of length 2 containing mean and standard deviation for
#' the log-normal distribution the alphas are generated from.
#' @param p_beta A vector of length 2 containing mean and standard deviation for
#' the normal distribution the alphas are generated from.
#' @param seed A seed for the random number generator.
#' @param params a logical value indicating wether true parameters should be
#' returned
#' 
#' @details
#' Generates data according to 
#' \deqn{
#'  \begin{equation}
#'  Y \sim Bernoulli(\pi_{jk}) \\
#'  pi_{jk} = (1 + exp(-alpha_k * (theta_j - beta_k)))^{-1}
#'  \end{equation}
#' }
#' 
#' @return If \code{params == TRUE} A list containing the generated data 
#' set and the true parameter vectors. 
#' 
#' If \code{params == FALSE} a matrix
#' containing the generated data
#' @export

sim_irt <- function(n_subjects, 
                    n_items,
                    p_theta = c(0, 1),
                    p_alpha = c(-1, 0.3),
                    p_beta = c(0, 1),
                    seed = NULL,
                    params = FALSE) {
  
  if(!is.null(seed)) set.seed(seed)

  theta <- rnorm(n_subjects, 0, 1)
  alpha <- rlnorm(n_items, -1, 0.35)
  beta  <- rnorm(n_items, 0, 1)
  tpar <- rbind(alpha, beta)
  
  # Draw data
  ip <- function(tpar, theta){
    1 / (1 + exp( - tpar[1] * (theta - tpar[2])))
  }
  # Probabilities
  ps <- apply(tpar, 2, ip, theta = theta)
  # Dichotomous data
  Y <- apply(ps, 2, rbinom, n = n_subjects, size = 1) 
  
  if(params) {
    return(list(Y = Y, theta = theta, alpha = alpha, beta = beta))
  } else {
    return(Y)
  }
}



