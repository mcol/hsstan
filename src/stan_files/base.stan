//
// Gaussian prior on regression coeffs (linear regression)
//

data {

  // number of unpenalized columns in model matrix
  int U;

  // number of training observations
  int N_train;

  // number of test observations
  int N_test;

  // prior standard deviation for the unpenalised variables
  int <lower=0> scale_u;

  // X matrix for training data
  matrix[N_train, U] X_train;

  // X matrix for test data
  matrix[N_test, U] X_test;

  // continuous response variable
  real y_train[N_train];
}

parameters {

  // unpenalized regression parameters
  vector[U] beta_u;

  // residual standard deviation
  real <lower=0> sigma;
}

model {

  // linear predictor
  vector[N_train] mu;
  mu = X_train[, 1:U] * beta_u;

  // unpenalized coefficients including intercept
  beta_u ~ normal(0, scale_u);

  // noninformative gamma priors on scale parameter are not advised
  sigma ~ inv_gamma(1, 1);

  // likelihood
  y_train ~ normal(mu, sigma);
}
