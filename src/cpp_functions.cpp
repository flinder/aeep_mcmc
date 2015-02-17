#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// [[Rcpp::depends("RcppArmadillo")]]

// multivariate normal density
// based on code by Ahmadou Dicko (http://gallery.rcpp.org/articles/dmvnorm_arma/)
// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma) {
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

// multivariate normal random number generator
// by Ahmadou Dicko (http://gallery.rcpp.org/articles/simulate-multivariate-normal/)
// [[Rcpp::export]]
arma::mat mvrnorm_arma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// Calculate the variance covariance matrix of the full posterior assuming
// multivariate normal sub-posteriors
// [[Rcpp::export(post_vcm)]]
arma::mat post_vcm(List post_list, int d, int M) {
  arma::mat vcm_post = arma::zeros(d, d);
  for(int i = 0; i < M; ++i) {
    arma::mat post = post_list[i];
    arma::mat vcm_m = arma::cov(post);
    arma::mat vcm_m_i = arma::inv(vcm_m);
    vcm_post += vcm_m_i;
  }
  arma::mat out = arma::inv(vcm_post);
  return(out);
}

// Calculate the mean vector of the full posterior assuming
// multivariate normal sub-posteriors
// [[Rcpp::export(post_mean)]]
arma::vec post_mean(List post_list, arma::mat post_vcm, int d, int M) {
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

// [[Rcpp::export(.theta_bar)]]
arma::vec theta_bar(arma::vec t, List post_list, double h, int d, int M) {
  arma::mat sel(M, d);
  for(int i = 0; i < M; ++i) {
    arma::mat postmat = post_list[i];
    sel.row(i) = postmat.row(t[i]);
  }
  arma::rowvec theta_bar = arma::mean(sel, 0);
  return trans(theta_bar);
}

// [[Rcpp::export(.mix_weight)]]
double mix_weight(arma::vec t, List post_list, double h, int d, int M) {
  arma::mat sel(M, d);
  for(int i = 0; i < M; ++i) {
    arma::mat postmat = post_list[i];
    sel.row(i) = postmat.row(t[i]);
  }
  arma::vec theta_b = theta_bar(t, post_list, h, d, M);
  arma::mat sigma = arma::eye(d, d) * pow(h, 2);
  arma::vec dens = dmvnrm_arma(sel, trans(theta_b), sigma);
  double out = prod(dens);
  return out;
}

//' Non-parametric Combination of Sub-posteriors
//' 
//' @param post_list A list containing the sub-posterior samples, stored in matrices
//' @return A matrix containing the samples from the sub-posterior product function
//' @export
// [[Rcpp::export]]
arma::mat combine_np(List post_list) {
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
  for(int i = 0; i < T; ++i) {
    double h = pow(i, (- 1 / (4 + d)));
    for(int m = 0; m < M; ++m) {
      c_dot = t_dot;
      c_dot[m] = urand[icount];
      double w_c_dot = mix_weight(c_dot, post_list, h, d, M);
      double w_t_dot = mix_weight(t_dot, post_list, h, d, M);
      double ratio = w_c_dot / w_t_dot;
      NumericVector u = runif(1);
      bool b = all(u < ratio).is_true();
      if(b) {
        t_dot = c_dot;
      }
      icount += 1;
    }
    arma::mat sigma = arma::eye(d, d) * pow(h, 2) / M;
    arma::vec mu = theta_bar(t_dot, post_list, h, d, M);
    arma::mat draws = mvrnorm_arma(1, mu, sigma);
    out.row(i) = draws;
  }
  return(out);
}

// Component weight function for semi parametric combination
// [[Rcpp::export(.mix_weight_sp)]]
arma::vec mix_weight_sp(arma::vec t, 
                     arma::mat sig_M, 
                     arma::vec mu_M, 
                     double w_t_dot,
                     List post_list,
                     double h,
                     int d,
                     int M,
                     List post_means,
                     List post_vcms) {
  arma::vec theta_b = theta_bar(t, post_list, h, d, M);
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

//' Semi-parametric Combination of Sub-posteriors
//' 
//' @param post_list A list containing the sub-posterior samples, stored in matrices
//' @return A matrix containing the samples from the sub-posterior product function
//' @export
// [[Rcpp::export]]
arma::mat combine_sp(List post_list) {
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
  arma::mat sig_M = post_vcm(post_list, d, M);
  arma::vec mu_M  = post_mean(post_list, sig_M, d, M);
  arma::vec sm_prod = sig_M.i() * mu_M;
  List post_means(M);
  List post_vcms(M);
  for(int m = 0; m < M; ++m) {
    arma::mat post = post_list[m];
    post_means[m] = mean(post);
    post_vcms[m] = arma::cov(post);
  }
  for(int i = 0; i < T; ++i) {
    double h = pow(i, (- 1 / (4 + d)));
    for(int m = 0; m < M; ++m) {
      c_dot = t_dot;
      c_dot[m] = urand[icount];
      double w_c_dot = mix_weight(c_dot, post_list, h, d, M);
      double w_t_dot = mix_weight(t_dot, post_list, h, d, M);
      arma::vec W_c_dot = mix_weight_sp(c_dot, sig_M, mu_M, w_c_dot, post_list, 
                                        h, d, M, post_means, post_vcms);
      arma::vec W_t_dot = mix_weight_sp(t_dot, sig_M, mu_M, w_t_dot, post_list, 
                                        h, d, M, post_means, post_vcms);
      arma::vec ratio_v = W_c_dot / W_t_dot;
      double ratio = ratio_v[0];
      NumericVector u_v = runif(1);
      double u = u_v[0];
      if(u < ratio) {
        t_dot = c_dot;
      }
      icount += 1;
    }
    arma::mat sig_t_dot = inv(arma::eye(d, d) * M / h + inv(sig_M));
    arma::vec theta_b = theta_bar(t_dot, post_list, h, d, M);
    arma::vec A = (M / h) * arma::eye(d, d) * theta_b + inv(sig_M) * mu_M;
    arma::vec mu_t_dot = sig_t_dot * A;
    arma::mat draws = mvrnorm_arma(1, mu_t_dot, sig_t_dot);
    out.row(i) = draws;
  }
  return(out);
}