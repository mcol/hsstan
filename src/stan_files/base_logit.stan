//
// Gaussian prior on regression coeffs (logistic regression)
//

data {

  // number of unpenalized columns in model matrix
  int U;

  // number of training observations
  int N_train;

  // number of test observations
  int N_test;

  // prior standard deviation for the unpenalised variables
  real <lower=0> scale_u;

  // X matrix for training data
  matrix[N_train, U] X_train;

  // X matrix for test data
  matrix[N_test, U] X_test;

  // binary response variable
  int<lower=0, upper=1> y_train[N_train];

  // binary response variable for test data
  int<lower=0, upper=1> y_test[N_test];
}

parameters {

  // unpenalized regression parameters
  vector[U] beta_u;
}

model {

  // linear predictor
  vector[N_train] mu = X_train * beta_u;

  // unpenalized coefficients including intercept
  beta_u ~ normal(0, scale_u);

  // likelihood
  y_train ~ bernoulli_logit(mu);
}
