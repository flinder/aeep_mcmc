#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]

// Rcpp::NumericMatrix to arma::mat
// [[Rcpp::export(.m_conv)]]
  arma::mat m_conv(NumericMatrix x) {
    arma::mat y = as<arma::mat>(x);
    return(y);
  }

// Rcpp::NumericVector to arma::rowvec
// [[Rcpp::export(.v_conv)]]
arma::mat v_conv(NumericVector x) {
  arma::mat y = as<arma::rowvec>(x);
  return(y);
}

// Rcpp::NumericVector to arma::vec
// [[Rcpp::export(.v_conv1)]]
arma::mat v_conv1(NumericVector x) {
  arma::mat y = as<arma::vec>(x);
  return(y);
}

// Mean from NumericVector
// [[Rcpp::export(.nv_mean)]]
double nv_mean(NumericVector x) {
  return sum(x) / x.size();
}

// multivariate normal density
// based on code by Ahmadou Dicko (http://gallery.rcpp.org/articles/dmvnorm_arma/)
// [[Rcpp::export]]
arma::vec dmvnrm_arma(NumericMatrix x_i,  
                      NumericVector mean_i,  
                      NumericMatrix sigma_i, 
                      bool logd = false) {
  const double log2pi = std::log(2.0 * M_PI);
  arma::mat x = m_conv(x_i);
  arma::mat sigma = m_conv(sigma_i);
  arma::rowvec mean = v_conv(mean_i);
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
  
  if (logd == false) {
    out = exp(out);
  }
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

// [[Rcpp::export(.theta_bar)]]
NumericVector theta_bar(IntegerVector t, List post_list, double h, int d, int M) {
  
  NumericMatrix sel(M, d);
  for(int i = 0; i < M; ++i) {
    NumericMatrix postmat = post_list[i];
    NumericMatrix::Row selrow = sel(i, _);
    NumericMatrix::Row pmrow = postmat(t[i], _);
    selrow = pmrow;
  }
  NumericVector theta_bar(d);
  for(int i = 0; i < sel.ncol(); ++i) {
    NumericMatrix::Column selcol = sel(_, i);
    theta_bar(i) = nv_mean(selcol);
  }
  return theta_bar;
}

// [[Rcpp::export(.mix_weight)]]
double mix_weight(IntegerVector t, List post_list, double h, int d, int M) {
  
  NumericMatrix sel(M, d);
  for(int i = 0; i < M; ++i) {
    NumericMatrix postmat = post_list[i];
    NumericMatrix::Row selrow = sel(i, _);
    NumericMatrix::Row pmrow = postmat(t[i], _);
    selrow = pmrow;
  }
  NumericVector theta_bar(d);
  for(int i = 0; i < sel.ncol(); ++i) {
    NumericMatrix::Column selcol = sel(_, i);
    theta_bar(i) = nv_mean(selcol);
  }
  NumericMatrix sigma(d, d);
  sigma.fill_diag(pow(h, 2));
  arma::vec dens = dmvnrm_arma(sel, theta_bar, sigma, FALSE);
  double out = prod(dens);
  return out;
}

// Calculate variance covariance matrix of a matrix
// [[Rcpp::export(.vcm)]]
arma::mat vcm(NumericMatrix X_i) {
  arma::mat X = m_conv(X_i);
  int n = X.n_rows;
  // transform inmat into deviation matrix
  arma::mat l = arma::ones(n, n);
  arma::mat x = X - l * X / n;
  arma::mat x_t = x.t();
  arma::mat V = (x_t * x) / n;
  // transform to Rcpp class
  //NumericMatrix out = wrap(V);
  return(V);
}

// [[Rcpp::export(post_vcm)]]
arma::mat post_vcm(List post_list) {
  NumericMatrix exmpl = post_list[0];
  int d = exmpl.ncol();
  int M = post_list.size();
  arma::mat vcm_post = arma::zeros(d, d);
  for(int i = 0; i < M; ++i) {
    arma::mat vcm_m = vcm(post_list[i]);
    arma::mat vcm_m_i = arma::inv(vcm_m);
    vcm_post += vcm_m_i;
  }
  arma::mat out = arma::inv(vcm_post);
  return(out);
}

// [[Rcpp::export(post_mean)]]
arma::vec post_mean(List post_list, arma::mat post_vcm) {
  NumericMatrix exmpl = post_list[0];
  //int n = exmpl.nrow();
  int d = exmpl.ncol();
  int M = post_list.size();
  arma::vec w_sig = arma::zeros(d, 1);
  for(int i = 0; i < M; ++i) {
    arma::mat vcm_m = vcm(post_list[i]);
    arma::mat vcm_m_i = arma::inv(vcm_m);
    arma::mat post_m = m_conv(post_list[i]);
    arma::rowvec mu_m = arma::mean(post_m, 0);
    arma::vec mu_m_t = mu_m.t();
    w_sig += vcm_m_i * mu_m_t;
  }
  arma::mat out = post_vcm * w_sig;
  //arma::vec out_t = out.t();
  return(out);
}

//' Non-parametric Combination of Sub-posteriors
//' 
//' @param post_list A list containing the sub-posterior samples, stored in matrices
//' @return A matrix containing the samples from the sub-posterior product function
// [[Rcpp::export]]
NumericMatrix combine_np(List post_list) {
  NumericMatrix exmpl = post_list[0];
  int d = exmpl.ncol();
  int T = exmpl.nrow();
  int M = post_list.size();
  IntegerVector Ts(T);
  for(int i = 0; i < T; i++) {
    Ts[i] = i;
  }
  IntegerVector t_dot = RcppArmadillo::sample(Ts, M, TRUE);  
  IntegerVector c_dot = t_dot;
  IntegerVector urand = RcppArmadillo::sample(Ts, (M * T), TRUE);
  int icount = 0;
  NumericMatrix out(T, d);
  
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
    NumericVector theta_b = theta_bar(t_dot, post_list, h, d, M);
    NumericMatrix sig(d, d);
    sig.fill_diag(pow(h, 2) / M);
    arma::vec mu = v_conv1(theta_b);
    arma::mat sigma = m_conv(sig);
    arma::mat draws = mvrnorm_arma(1, mu, sigma);
    NumericVector draws_rcpp = wrap(draws);
    NumericMatrix::Row outrow = out(i, _);
    outrow = draws_rcpp;
	}
  return(out);
}


//' Semi-parametric Combination of Sub-posteriors
//' 
//' @param post_list A list containing the sub-posterior samples, stored in matrices
//' @return A matrix containing the samples from the sub-posterior product function
// [[Rcpp::export]]
NumericMatrix combine_sp(List post_list) {
  NumericMatrix exmpl = post_list[0];
  int d = exmpl.ncol();
  int T = exmpl.nrow();
  int M = post_list.size();
  IntegerVector Ts(T);
  for(int i = 0; i < T; i++) {
    Ts[i] = i;
  }
  IntegerVector t_dot = RcppArmadillo::sample(Ts, M, TRUE);  
  IntegerVector c_dot = t_dot;
  IntegerVector urand = RcppArmadillo::sample(Ts, (M * T), TRUE);
  int icount = 0;
  NumericMatrix out(T, d);
  arma::mat sig_m = post_vcm(post_list);
  arma::vec mu_m  = post_mean(post_list, sig_m);
  arma::vec sm_prod = sig_m.i() * mu_m;
  
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
    NumericVector theta_b = theta_bar(t_dot, post_list, h, d, M);
    arma::vec theta_b_con = v_conv1(theta_b);
    arma::mat sig(d, d);
    sig.eye();
    sig *= (M / h);
    arma::mat sig1(d, d);
    sig1 += sig_m.i();
    arma::mat sig_t_dot = sig1.i();
    arma::vec mu_t_dot = sig_t_dot * (sig * theta_b_con + sm_prod);
    arma::mat draws = mvrnorm_arma(1, theta_b_con, sig_t_dot);
    NumericVector draws_rcpp = wrap(draws);
    NumericMatrix::Row outrow = out(i, _);
    outrow = draws_rcpp;
	}
  return(out);
}