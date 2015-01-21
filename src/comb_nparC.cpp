#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]

// Multivariate Normal Density
const double log2pi = std::log(2.0 * M_PI);

// Rcpp::NumericMatrix to arma::mat
  arma::mat m_conv(NumericMatrix x) {
    arma::mat y = as<arma::mat>(x);
    return(y);
  }

// Rcpp::NumericVector to arma::rowvec
arma::mat v_conv(NumericVector x) {
  arma::mat y = as<arma::rowvec>(x);
  return(y);
}

arma::mat v_conv1(NumericVector x) {
  arma::mat y = as<arma::vec>(x);
  return(y);
}

// Mean from NumericVector
double nv_mean(NumericVector x) {
  return sum(x) / x.size();
}

// multivariate normal density
// by Ahmadou Dicko (http://gallery.rcpp.org/articles/dmvnorm_arma/)
arma::vec dmvnrm_arma(NumericMatrix x_i,  
                      NumericVector mean_i,  
                      NumericMatrix sigma_i, 
                      bool logd = false) {
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
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


NumericVector theta_bar(IntegerVector t, List mcmcout, double h, int d, int M) {
  NumericMatrix sel(M, d);
  for(int i = 0; i < M; ++i) {
    NumericMatrix postmat = mcmcout[i];
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

double mix_weight(IntegerVector t, List mcmcout, double h, int d, int M) {
  NumericMatrix sel(M, d);
  for(int i = 0; i < M; ++i) {
    NumericMatrix postmat = mcmcout[i];
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

/* Main Function */

// [[Rcpp::export]]
NumericMatrix comb_nparC(List mcmcout) {
  
  NumericMatrix expl = mcmcout[0];
  int d = expl.ncol();
  int T = expl.nrow();
	int M = mcmcout.size();
  //Rcpp::Rcout << "d = " << d << std::endl;
  //Rcpp::Rcout << "T = " << T << std::endl;
  //Rcpp::Rcout << "M = " << M << std::endl;
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
      double w_c_dot = mix_weight(c_dot, mcmcout, h, d, M);
      double w_t_dot = mix_weight(t_dot, mcmcout, h, d, M);
      double ratio = w_c_dot / w_t_dot;
      NumericVector u = runif(1);
      bool b = all(u < ratio).is_true();
      if(b) {
        t_dot = c_dot;
      }
      icount += 1;
		}
    NumericVector theta_b = theta_bar(t_dot, mcmcout, h, d, M);
    NumericMatrix sig(d, d);
    sig.fill_diag(pow(h, 2) / M);
    arma::vec mu = v_conv1(theta_b);
    arma::mat sigma = m_conv(sig);
    arma::mat draws = mvrnormArma(1, mu, sigma);
    NumericVector draws_rcpp = wrap(draws);
    NumericMatrix::Row outrow = out(i, _);
    outrow = draws_rcpp;
	}
  return(out);
}