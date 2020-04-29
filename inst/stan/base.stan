//
// Gaussian prior on regression coeffs (linear regression)
//

data {

  // number of unpenalized columns in model matrix
  int U;

  // number of observations
  int N;

  // design matrix
  matrix[N, U] X;

  // continuous response variable
  vector[N] y;

  // prior standard deviation for the unpenalised variables
  real<lower=0> scale_u;
}

parameters {

  // unpenalized regression parameters
  vector[U] beta_u;

  // residual standard deviation
  real <lower=0> sigma;
}

model {

  // unpenalized coefficients including intercept
  beta_u ~ normal(0, scale_u);

  // noninformative gamma priors on scale parameter are not advised
  sigma ~ inv_gamma(1, 1);

  // likelihood
  y ~ normal_id_glm(X, 0, beta_u, sigma);
}
