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

transformed parameters {

  // linear predictor
  vector[N_train] theta_train;

  theta_train = X_train[, 1:U] * beta_u;
}

model {

  // unpenalized coefficients including intercept
  beta_u ~ normal(0, 1000);

  // noninformative gamma priors on scale parameter are not advised
  sigma ~ inv_gamma(1, 1);

  // likelihood
  y_train ~ normal(theta_train, sigma);
}

generated quantities {

  // predicted outcome on widthdrawn data
  vector[N_test] y_pred = X_test[, 1:U] * beta_u;
}
