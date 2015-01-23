functions {
  real dw_normal_log(real y, real mu, real sigma, int M) {
    return normal_log(y, mu, sigma) / M;
  }
  real dw_lognormal_log(real y, real mu, real sigma, int M) {
    return lognormal_log(y, mu, sigma) / M;
  }
}
data {
  int<lower=1> J; // number of students
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of observations
  int<lower=1,upper=J> jj[N]; // student for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> Y[N]; // correctness for observation n
  int<lower=1> M; // number of partitions
}
parameters {
  real<lower=0, upper=5> alpha[K]; // discrimination parameters
  real beta[K]; // difficulty parameters
  real theta[J]; // student ability parameters
}
model {
  real pi[N];
  for(j in 1:J)
    theta[j] ~ dw_normal(0.0, 1.0, M);
  for(k in 1:K) {
    beta[k] ~ dw_normal(0.0, 1.0, M);
    alpha[k] ~ dw_lognormal(0.0, 6.0, M);    
  }
 
  for (n in 1:N)
    pi[n] <- 1 / (1 + exp( - alpha[kk[n]] * (theta[jj[n]] - beta[kk[n]])));
  
  Y ~ bernoulli(pi);
}