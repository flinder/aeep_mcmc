#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// [[Rcpp::depends("RcppArmadillo")]]

// Where appropriate I reference to Neiswanger et al. 2014 
// e.g. eq 3.1 refers to equation 3.1 in the paper
//		p 3 refers to page 3

// multivariate normal density
// Returns a scalar density for each row of X
// given mean vector 'mean' and variance covariance matrix 'sigma'
// based on code by Ahmadou Dicko (http://gallery.rcpp.org/articles/dmvnorm_arma/)
// [[Rcpp::export]]
arma::vec dmvnrm_arma(const arma::mat& x,  
                      const arma::rowvec& mean,  
                      const arma::mat& sigma) {
  const double log2pi = std::log(2.0 * M_PI);
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
 
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  out = exp(out);
  return(out);
}

// Multivariate normal random number generator
// Returns n random vectors of size mu from the multivariate normal 
// with mean mu and variance covariance matrix sigma
// by Ahmadou Dicko (http://gallery.rcpp.org/articles/simulate-multivariate-normal/)
// [[Rcpp::export]]
arma::mat mvrnorm_arma(int n, const arma::vec& mu, const arma::mat& sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// Calculate the variance covariance matrix of the full posterior assuming
// multivariate normal sub-posteriors
// eq 3.1
// Input: 
// post_list: List of M sub-posteriors (n*d matrices)
// M: Number of sub posteriors
// d: dimensionality of the posteriors
// [[Rcpp::export(post_vcm)]]
arma::mat post_vcm(const List& post_list, int d, int M) {
  arma::mat vcm_post = arma::zeros(d, d);
  for(int i = 0; i < M; ++i) {
    arma::mat post = post_list[i];
    arma::mat vcm_m = arma::cov(post); // calculate variance covariance matrix
    arma::mat vcm_m_i = arma::inv(vcm_m); // invert it
    vcm_post += vcm_m_i;
  }
  arma::mat out = arma::inv(vcm_post);
  return(out);
}

// Calculate the mean vector of the full posterior assuming
// multivariate normal sub-posteriors
// eq. 3.2
// Input:
// post_list: List of M sub-posteriors (n*d matrices)
// M: Number of sub posteriors
// d: dimensionality of the posteriors
// post_vcm: combined variance covariance matrix (output of post_vcm())
// [[Rcpp::export(post_mean)]]
arma::vec post_mean(const List& post_list, const arma::mat& post_vcm, int d, int M) {
  arma::vec w_sig = arma::zeros(d, 1);
  for(int i = 0; i < M; ++i) {
    arma::mat post = post_list[i];
    arma::mat vcm_m = arma::cov(post);
    arma::rowvec mu_m = mean(post, 0);
    w_sig += inv(vcm_m) * trans(mu_m);
  }
  arma::mat out = post_vcm * w_sig;
  return(out);
}

// Calculate the mean vector for selected sample in iteration of sampler for
// non- and semi-parametric combination
// eq 3.4
// Input:
// t: indices of currently selected samples
// post_list: List of sub-posteriors
// h: Bandwidth parameter in the current state
// d: dimensionality of the posterior (parameters in the model)
// M: number of sub-posteriors
// [[Rcpp::export(.theta_bar)]]
arma::vec theta_bar(const arma::vec& t, const List& post_list, double h, int d, int M) {
  // Create matrix of selected samples
  arma::mat sel(M, d);
  for(int i = 0; i < M; ++i) {
    arma::mat postmat = post_list[i];
    sel.row(i) = postmat.row(t[i]);
  }
  arma::rowvec theta_bar = arma::mean(sel, 0);
  return trans(theta_bar);
}

// Calculates mixture weights of current component 
// for semi and non parametric combination
// eq 3.5
// Input:
// t: indices of currently selected samples
// post_list: List of sub-posteriors
// h: Bandwidth parameter in the current state
// d: dimensionality of the posterior (parameters in the model)
// M: number of sub-posteriors
// [[Rcpp::export(.mix_weight)]]
double mix_weight(const arma::vec& t, 
                  const List& post_list, 
                  double h, 
                  int d, 
                  int M, 
                  const arma::vec& theta_b) {
  
  // Create matrix of selected samples
  arma::mat sel(M, d);
  for(int i = 0; i < M; ++i) {
    arma::mat postmat = post_list[i];
    sel.row(i) = postmat.row(t[i]);
  }
  arma::mat sigma = arma::eye(d, d) * pow(h, 2);
  // Calculate density of selected samples
  arma::vec dens = dmvnrm_arma(sel, trans(theta_b), sigma);
  double out = prod(dens);
  return out;
}

// Algorithm 1 on page 5

// This stuff is for the R Package documentation (Roxygen2):

//' Non-parametric Combination of Sub-posteriors
//' @param post_list A list containing the sub-posterior samples, stored in matrices
//' @return A matrix containing the samples from the sub-posterior product function
//' @export
// [[Rcpp::export]]
arma::mat combine_np(const List& post_list) {
  // Get number of dimensions , sub-posteriors, etc...
  arma::mat exmpl = post_list[0];
  int d = exmpl.n_cols; // dimensions
  int T = exmpl.n_rows; // # of mcmc iterations per sub-posterior
  int M = post_list.size(); // # of sub-posteriors
  arma::vec Ts(T);
  // This just creates a sequence from 1 to T
  for(int i = 0; i < T; i++) {
    Ts[i] = i;
  }
  // sample first set of indices
  arma::vec t_dot = RcppArmadillo::sample(Ts, M, TRUE);  
  arma::vec c_dot = t_dot;
  // sample uniform random numbers to generate new selection indices
  arma::vec urand = RcppArmadillo::sample(Ts, (M * T), TRUE);
  int icount = 0; //index to loop through urand
  arma::mat out(T, d); // output matrix
  
  // Gibbs loop
  for(int i = 0; i < T; ++i) {
    double h = pow(i, (- 1 / (4 + d))); // bandwidth parameter
    
    // Metropolis Loop
    for(int m = 0; m < M; ++m) { 
      
      c_dot = t_dot;
      // change one index for proposal
      c_dot[m] = urand[icount];
      
      // calculate mixture weight of proposal
      arma::vec theta_b = theta_bar(c_dot, post_list, h, d, M);
      double w_c_dot = mix_weight(c_dot, post_list, h, d, M, theta_b);
      
      // calculate old mixture weight (this step migt be avoidable by carrying 
      // the weight calculated above to the next iteration)
      theta_b = theta_bar(t_dot, post_list, h, d, M);
      double w_t_dot = mix_weight(t_dot, post_list, h, d, M, theta_b);
      
      // get density ratio of proposal / last iteration
      double ratio = w_c_dot / w_t_dot;
      NumericVector u = runif(1);
      
      // Metropolis stochastic acceptance step
      bool b = all(u < ratio).is_true();
      if(b) {
        t_dot = c_dot;
      }
      icount += 1;
    }
    // Draw from mixture in gibbs step
    arma::mat sigma = arma::eye(d, d) * pow(h, 2) / M;
    arma::vec mu = theta_bar(t_dot, post_list, h, d, M);
    arma::mat draws = mvrnorm_arma(1, mu, sigma);
    out.row(i) = draws;
  }
  return(out);
}

// Component weight function for semi parametric combination
// as described on p 6 (W_t_dot)
// Input:
// t: current selection
// sig_M: Parametric approximation to variance-covariance matrix of posterior
//  (see function post_vcm() above)
// mu_M: Parametric approxiation to posterior mean (see function post_mean() above)
// w_t_dot: component weight (see function mix_weight())
// h, d, M: as in functions above
// post_means: list of mean vectors of sub-posteriors
// post_vcms: Variance covariance matrices of sub-posteriors
// theta_b: Mean vector in current iteration
// [[Rcpp::export(.mix_weight_sp)]]
arma::vec mix_weight_sp(const arma::vec& t, 
                     const arma::mat& sig_M, 
                     const arma::vec& mu_M, 
                     const double& w_t_dot,
                     const List& post_list,
                     const double& h,
                     const int& d,
                     const int& M,
                     const List& post_means,
                     const List& post_vcms,
                     const arma::vec& theta_b) {
  arma::mat theta_b_mat(1, d);
  theta_b_mat.row(0) = trans(theta_b);
  arma::mat sig = sig_M + arma::eye(d, d) * h / M;
  arma::vec num = w_t_dot * dmvnrm_arma(theta_b_mat, trans(mu_M), sig);
  arma::vec den = arma::ones(1);
  for(int m = 0; m < M; ++m) {
    arma::mat post = post_list[m];
    arma::rowvec theta_tm = post.row(t[m]);
    arma::rowvec mu_m = post_means[m];
    arma::mat vcm_m = post_vcms[m];
    arma::vec d = dmvnrm_arma(theta_tm, mu_m, vcm_m);
    den = den * d;
  }
  arma::vec out = num / den;
  return out;
}

// Semi parametric combination as described in section 3.3. Very similar to 
// non-parametric combination. 

// This stuff is for R package documentation:
//' Semi-parametric Combination of Sub-posteriors
//' 
//' @param post_list A list containing the sub-posterior samples, stored in matrices
//' @return A matrix containing the samples from the sub-posterior product function
//' @export
// [[Rcpp::export]]
arma::mat combine_sp(const List& post_list) {
  arma::mat exmpl = post_list[0];
  int d = exmpl.n_cols;
  int T = exmpl.n_rows;
  int M = post_list.size();
  arma::vec Ts(T);
  for(int i = 0; i < T; i++) {
    Ts[i] = i;
  }
  arma::vec t_dot = RcppArmadillo::sample(Ts, M, TRUE);  
  arma::vec c_dot = t_dot;
  arma::vec urand = RcppArmadillo::sample(Ts, (M * T), TRUE);
  int icount = 0;
  arma::mat out(T, d);
  
  // Calculate parametric approximations of full posterior and vcm
  arma::mat sig_M = post_vcm(post_list, d, M);
  arma::vec mu_M  = post_mean(post_list, sig_M, d, M);
  arma::vec sm_prod = sig_M.i() * mu_M;
  
  // Calculate means and vcms for sub-posteriors to avoid doing it in each iteration
  List post_means(M);
  List post_vcms(M);
  for(int m = 0; m < M; ++m) {
    arma::mat post = post_list[m];
    post_means[m] = mean(post);
    post_vcms[m] = arma::cov(post);
  }
  
  // Gibbs Loop
  for(int i = 0; i < T; ++i) {
    double h = pow(i, (- 1 / (4 + d)));
    
    // Metropolis Loop
    for(int m = 0; m < M; ++m) {
      c_dot = t_dot;
      // Proposal
      c_dot[m] = urand[icount];
      
      // Mixture weight for proposal
      arma::vec theta_b = theta_bar(c_dot, post_list, h, d, M);
      double w_c_dot = mix_weight(c_dot, post_list, h, d, M, theta_b);
      arma::vec W_c_dot = mix_weight_sp(c_dot, sig_M, mu_M, w_c_dot, post_list, 
                                        h, d, M, post_means, post_vcms, theta_b);
      
      // Mixture weight for last step (again can be avoided by using above calculated 
      // weight from the last iteration?)
      theta_b = theta_bar(t_dot, post_list, h, d, M);
      double w_t_dot = mix_weight(t_dot, post_list, h, d, M, theta_b);
      arma::vec W_t_dot = mix_weight_sp(t_dot, sig_M, mu_M, w_t_dot, post_list, 
                                        h, d, M, post_means, post_vcms, theta_b);
      // Density ratio
      arma::vec ratio_v = W_c_dot / W_t_dot;
      
      // acceptance step
      double ratio = ratio_v[0];
      NumericVector u_v = runif(1);
      double u = u_v[0];
      if(u < ratio) {
        t_dot = c_dot;
      }
      icount += 1;
    }
    
    // Draw Gibbs samples
    arma::mat sig_t_dot = inv(arma::eye(d, d) * M / h + inv(sig_M));
    arma::vec theta_b = theta_bar(t_dot, post_list, h, d, M);
    arma::vec A = (M / h) * arma::eye(d, d) * theta_b + inv(sig_M) * mu_M;
    arma::vec mu_t_dot = sig_t_dot * A;
    arma::mat draws = mvrnorm_arma(1, mu_t_dot, sig_t_dot);
    out.row(i) = draws;
  }
  return(out);
}
