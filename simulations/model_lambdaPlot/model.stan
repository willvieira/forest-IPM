data {
  int<lower=1> N;
  vector [N] lambda_log;
  vector [N] MAT;
}
parameters {
  // slopes
  real beta_mean; 
  real beta_sigma;
  // intercepts
  real inter;
  real sigma_inter;
  // shape skewness
  real alpha;
}
model {
  // Likelihood
  lambda_log ~ skew_normal(
    beta_mean * MAT + inter,                // mean
    exp(sigma_inter + beta_sigma * MAT),    // SD
    alpha                                   // shape (skewness intensity)
  );
  
  // Priors
  beta_mean ~ normal(0, 1);
  beta_sigma ~ normal(0, 1);
  inter ~ normal(0, 1);
  sigma_inter ~ normal(0, 1);
  alpha ~ normal(0, 1);
}
