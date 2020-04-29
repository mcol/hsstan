//
// Gaussian prior on regression coeffs (logistic regression)
//

data {

  // number of unpenalized columns in model matrix
  int U;

  // number of observations
  int N;

  // prior standard deviation for the unpenalised variables
  real <lower=0> scale_u;

  // design matrix
  matrix[N, U] X;

  // binary response variable
  int<lower=0, upper=1> y[N];
}

parameters {

  // unpenalized regression parameters
  vector[U] beta_u;
}

model {

  // unpenalized coefficients including intercept
  beta_u ~ normal(0, scale_u);

  // likelihood
  y ~ bernoulli_logit_glm(X, 0, beta_u);
}
