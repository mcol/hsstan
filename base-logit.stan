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

  // X matrix for training data
  matrix[N_train, U] X_train;

  // X matrix for test data
  matrix[N_test, U] X_test;

  // binary response variable
  int<lower=0, upper=1> y_train[N_train];
}

parameters {

  // unpenalized regression parameters
  vector[U] beta_u;
}

transformed parameters {

  // linear predictor
  vector[N_train] theta_train;

  theta_train = X_train[, 1:U] * beta_u;
}

model {

  // unpenalized coefficients including intercept
  beta_u ~ normal(0, 1000);

  y_train ~ bernoulli_logit(theta_train);
}

generated quantities {

  // predicted outcome on widthdrawn data
  vector[N_test] y_pred = X_test[, 1:U] * beta_u;
}
